# Copyright (C) 2018 Alejandro Jaramillo
# Copyright (C) 2026 Alejandro Jaramillo-Moreno
# SPDX-License-Identifier: GPL-3.0-only
"""
PyWKSpectra - Wheeler-Kiladis Spectral Analysis
===============================================

Wheeler-Kiladis Space-Time Spectral Analysis in Python

Author: Alejandro Jaramillo-Moreno
        Departamento de Cambio Climático
        Instituto de Ciencias de la Atmósfera y Cambio Climático
        Universidad Nacional Autónoma de México (UNAM)

Description
-----------
Scientific functions for computing Wheeler-Kiladis space-time
spectra, including temporal preprocessing, equatorial symmetry
decomposition, Fourier analysis, spectral power estimation,
and background normalization.
"""

import numpy as np
import xarray as xr

from .io import load_wk_data

# ==================================================
# 1. Temporal sampling parameters
# ==================================================
def sampling_vars(array_in, spd, nDayWin, nDaySkip, verbose=True):
    """
    Calculate sampling parameters for Wheeler-Kiladis spectral analysis.

    Parameters
    ----------
    array_in : numpy.ndarray or xarray.DataArray
        Input data with dimensions (time, lat, lon).
    spd : int
        Number of samples per day.
    nDayWin : int
        Spectral window length in days.
    nDaySkip : int
        Number of days skipped between consecutive windows.
        Negative values indicate overlap.
        Zero indicates consecutive, non-overlapping windows.
        Positive values indicate gaps between windows.
    verbose : bool, default True
        Print a summary of the sampling configuration.

    Returns
    -------
    ntim, nlat, nlon : int
        Input array dimensions.
    nDayTot : int
        Total number of complete days.
    nSampWin : int
        Number of samples per spectral window.
    nSampSkip : int
        Number of skipped samples between windows.
    """

    # Validate input dimensions
    if array_in.ndim != 3:
        raise ValueError(
            "Input must have dimensions (time, lat, lon)."
        )

    # Validate sampling parameters
    for name, value in [
        ("spd", spd),
        ("nDayWin", nDayWin),
        ("nDaySkip", nDaySkip)
    ]:
        if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, np.integer)
        ):
            raise ValueError(f"{name} must be an integer.")

    if spd <= 0 or nDayWin <= 0:
        raise ValueError(
            "spd and nDayWin must be positive integers."
        )

    # Input dimensions
    ntim, nlat, nlon = array_in.shape

    # Require complete days
    if ntim % spd != 0:
        raise ValueError(
            "Input must contain complete days. "
            f"Remaining samples: {ntim % spd}"
        )

    # Sampling parameters
    nDayTot = ntim // spd
    nSampWin = nDayWin * spd
    nSampSkip = nDaySkip * spd

    # Validate window configuration
    if nSampWin > ntim:
        raise ValueError(
            "Spectral window exceeds the available time series."
        )

    nSampStep = nSampWin + nSampSkip

    if nSampStep <= 0:
        raise ValueError(
            "Window advance must be greater than zero."
        )

    # Number of complete spectral windows
    nWindows = 1 + (ntim - nSampWin) // nSampStep

    # Unused samples at the end
    nUnused = (ntim - nSampWin) % nSampStep

    # Print information
    if verbose:

        print("\nPyWKSpectra - Spectral Sampling")
        print("=" * 55)

        print(f"Time samples:       {ntim}")
        print(f"Latitude points:    {nlat}")
        print(f"Longitude points:   {nlon}")
        print(f"Samples per day:    {spd}")
        print(f"Total days:         {nDayTot}")

        print("-" * 55)

        print(f"Window length:      {nDayWin} days")
        print(f"Window samples:     {nSampWin}")
        print(f"Window skip:        {nDaySkip} days")
        print(f"Window advance:     {nSampStep / spd:g} days")
        print(f"Window overlap:     {max(0, -nDaySkip)} days")

        print("-" * 55)

        print(f"Complete windows:   {nWindows}")
        print(f"Unused final days:  {nUnused / spd:g}")

        print("=" * 55)

    return ntim, nlat, nlon, nDayTot, nSampWin, nSampSkip

# ==================================================
# 2. Remove dominant low-frequency signals
# ==================================================
def remove_dominant_signals(array_in, spd, nDayTot, nSampWin, progress=True, lat_chunk=2):
    """
    Remove dominant low-frequency variability from atmospheric data.

    A linear trend is removed while preserving the temporal mean.
    For records of at least 365 days, Fourier components with
    frequencies 0 < f <= 1/window_days are also removed.

    Parameters
    ----------
    array_in : numpy.ndarray or xarray.DataArray
        Input data with dimensions (time, lat, lon).
    spd : int
        Number of samples per day.
    nDayTot : int
        Total number of days.
    nSampWin : int
        Number of samples per spectral window.
    progress : bool, default True
        Display processing information and progress bar.
    lat_chunk : int, default 2
        Number of latitudes processed simultaneously.

    Returns
    -------
    array_dt : numpy.ndarray or xarray.DataArray
        Processed data with the same dimensions and coordinates.

    Notes
    -----
    The temporal mean is preserved. Low-frequency removal uses
    a sharp Fourier cutoff over the complete time series.
    """

    # Input dimensions
    ntim, nlat, nlon = array_in.shape

    # Centered time coordinate
    t = np.arange(ntim, dtype=np.float64)
    t -= t.mean()
    t2 = np.dot(t, t)

    # Critical frequency (cycles/day)
    fCrit = spd / nSampWin
    apply_filter = nDayTot >= 365

    # Fourier frequencies
    if apply_filter:
        freq = np.fft.rfftfreq(ntim, d=1.0 / spd)
        remove = (freq > 0) & (freq <= fCrit)

    # Initialize output
    output_dtype = np.float32 if array_in.dtype == np.float32 else np.float64
    result = np.empty((ntim, nlat, nlon), dtype=output_dtype)

    # Processing information
    if progress:
        print("\nPyWKSpectra - Removing Dominant Signals")
        print("=" * 55)
        print("Linear trend:       Removed")
        print("Temporal mean:      Preserved")
        print(f"High-pass filter:   {'Applied' if apply_filter else 'Skipped'}")

        if apply_filter:
            print(f"Critical frequency: {fCrit:.5f} cycles/day")
            print(f"Removed periods:    >= {nSampWin / spd:g} days")

        print(f"Latitude chunk:     {lat_chunk}")
        print("-" * 55)

    # Process latitude blocks
    for start in range(0, nlat, lat_chunk):

        stop = min(start + lat_chunk, nlat)

        # Extract block
        block = np.array(array_in[:, start:stop, :], dtype=np.float64, copy=True)

        # Remove linear trend while preserving temporal mean
        values = block.reshape(ntim, -1)
        slope = (t @ values) / t2
        values -= t[:, None] * slope[None, :]

        # Remove low-frequency Fourier components
        if apply_filter:
            fourier = np.fft.rfft(block, axis=0)
            fourier[remove, :, :] = 0.0
            block = np.fft.irfft(fourier, n=ntim, axis=0)

        # Store processed data
        result[:, start:stop, :] = block

        # Progress bar
        if progress:
            percentage = 100 * stop / nlat
            filled = round(30 * stop / nlat)
            bar = "#" * filled + "-" * (30 - filled)
            print(f"\rProgress: [{bar}] {percentage:5.1f}%", end="", flush=True)

    if progress:
        print("\nProcessing completed.")
        print("=" * 55)

    # Preserve xarray structure
    if isinstance(array_in, xr.DataArray):
        array_dt = xr.DataArray(result, dims=array_in.dims, coords=array_in.coords, attrs=array_in.attrs.copy(), name=array_in.name)
    else:
        array_dt = result

    return array_dt

# ==================================================
# 3. Equatorial symmetry decomposition
# ==================================================
def decompose_symasym(array_in, verbose=True):
    """
    Decompose atmospheric data into equatorially symmetric and
    antisymmetric components.

    Symmetric components are stored in the Southern Hemisphere,
    and antisymmetric components in the Northern Hemisphere.
    Equatorial values remain unchanged.

    Parameters
    ----------
    array_in : xarray.DataArray
        Input data with dimensions (time, lat, lon).
        Latitudes must be ordered north to south and symmetric
        about the equator, as provided by load_wk_data().
    verbose : bool, default True
        Print decomposition information.

    Returns
    -------
    array_as : xarray.DataArray
        Combined symmetric and antisymmetric components with
        the same dimensions and coordinates as the input.
    """

    # Mirror data across the equator
    mirrored = array_in.isel(lat=slice(None, None, -1)).assign_coords(lat=array_in.lat)

    # Symmetric (south) and antisymmetric (north) components
    array_as = xr.where(array_in.lat > 0, 0.5 * (array_in - mirrored), 0.5 * (array_in + mirrored))

    # Preserve metadata
    array_as = array_as.transpose(*array_in.dims)
    array_as.name = array_in.name
    array_as.attrs = array_in.attrs.copy()

    # Processing information
    if verbose:
        print("\nPyWKSpectra - Equatorial Decomposition")
        print("=" * 55)
        print("Symmetric:       Southern Hemisphere")
        print("Antisymmetric:   Northern Hemisphere")
        print(f"Latitude pairs:  {np.count_nonzero(array_in.lat.values > 0)}")
        print(f"Equator present: {np.any(np.isclose(array_in.lat.values, 0))}")
        print(f"Output shape:    {array_as.shape}")
        print("=" * 55)

    return array_as

# ==================================================
# 4. Space-time Fourier coefficients
# ==================================================
def spectral_coefficients(array_in, spd, nSampWin, nSampSkip, progress=True):
    """
    Calculate space-time Fourier coefficients for Wheeler-Kiladis analysis.

    Each temporal window is linearly detrended, multiplied by a Hann
    taper, and transformed using a two-dimensional FFT along time
    and longitude.

    Fourier coefficients are reordered following the Hayashi
    convention used in the original wk_spectra implementation.

    Parameters
    ----------
    array_in : numpy.ndarray or xarray.DataArray
        Input data with dimensions (time, lat, lon).
    spd : int
        Number of samples per day.
    nSampWin : int
        Number of samples per spectral window.
    nSampSkip : int
        Samples skipped between windows. Negative values indicate overlap.
    progress : bool, default True
        Display processing information and progress bar.

    Returns
    -------
    wavefft : numpy.ndarray
        Zonal wavenumbers.
    freqfft : numpy.ndarray
        Temporal frequencies in cycles per day.
    peeAS : numpy.ndarray
        Complex Fourier coefficients with dimensions
        (window, frequency, latitude, wavenumber).

    Notes
    -----
    - Uses a Hann temporal window.
    - FFT coefficients are normalized by nSampWin * nlon.
    - Nyquist endpoints are duplicated for even dimensions
      to preserve compatibility with the original implementation.
    - No spectral power or latitude averaging is calculated.
    """

    # Input dimensions
    ntim, nlat, nlon = array_in.shape

    # Temporal window configuration
    nSampStep = nSampWin + nSampSkip

    if nSampWin < 3 or nSampWin > ntim or nSampStep <= 0:
        raise ValueError("Invalid spectral window configuration.")

    # Number of complete windows
    nWindow = (ntim - nSampWin) // nSampStep + 1

    # Fourier frequency and wavenumber coordinates
    freqfft = np.fft.fftshift(np.fft.fftfreq(nSampWin, d=1.0 / spd))
    wavefft = np.arange(-(nlon // 2), nlon // 2 + 1)

    # Include positive Nyquist endpoint for even time dimension
    if nSampWin % 2 == 0:
        freqfft = np.append(freqfft, spd / 2)

    # Allocate complex Fourier coefficients
    peeAS = np.empty((nWindow, len(freqfft), nlat, len(wavefft)), dtype=np.complex128)

    # Centered time coordinate for linear detrending
    t = np.arange(nSampWin, dtype=np.float64)
    t -= t.mean()
    t2 = np.dot(t, t)

    # Hann temporal taper
    taper = np.hanning(nSampWin)[:, None, None]

    # FFT normalization
    norm = nSampWin * nlon

    # Processing information
    if progress:
        print("\nPyWKSpectra - Spectral Coefficients")
        print("=" * 55)
        print(f"Spectral windows:   {nWindow}")
        print(f"Window length:      {nSampWin / spd:g} days")
        print(f"Window overlap:     {max(0, -nSampSkip) / spd:g} days")
        print(f"Window advance:     {nSampStep / spd:g} days")
        print("Linear detrending:  Applied")
        print("Temporal taper:     Hann")
        print("FFT dimensions:     Time and longitude")
        print(f"Output memory:      {peeAS.nbytes / 1024**2:.2f} MiB")
        print("-" * 55)

    # Process spectral windows
    for nw in range(nWindow):

        # Temporal window indices
        start = nw * nSampStep
        stop = start + nSampWin

        # Extract temporal window
        window = np.array(array_in[start:stop, :, :], dtype=np.float64, copy=True)

        # Remove temporal mean and linear trend
        values = window.reshape(nSampWin, -1)
        values -= values.mean(axis=0, keepdims=True)
        slope = (t @ values) / t2
        values -= t[:, None] * slope[None, :]

        # Apply Hann temporal taper
        window *= taper

        # Space-time Fourier transform
        fourier = np.fft.fft2(window, axes=(0, 2)) / norm

        # Center temporal frequencies and zonal wavenumbers
        fourier = np.fft.fftshift(fourier, axes=(0, 2))

        # Reorder zonal wavenumbers (Hayashi convention)
        if nlon % 2 == 0:
            peeAS[nw, :nSampWin, :, 0] = fourier[:, :, 0]
            peeAS[nw, :nSampWin, :, 1:] = fourier[:, :, ::-1]
        else:
            peeAS[nw, :nSampWin, :, :] = fourier[:, :, ::-1]

        # Duplicate temporal Nyquist endpoint when required
        if nSampWin % 2 == 0:
            peeAS[nw, nSampWin, :, :] = peeAS[nw, 0, :, :]

        # Update progress bar
        if progress:
            percentage = 100 * (nw + 1) / nWindow
            filled = round(30 * (nw + 1) / nWindow)
            bar = "#" * filled + "-" * (30 - filled)
            print(f"\rProgress: [{bar}] {percentage:5.1f}% ({nw + 1}/{nWindow})", end="", flush=True)

    if progress:
        print("\nProcessing completed.")
        print("=" * 55)

    return wavefft, freqfft, peeAS

# ==================================================
# 5. Spectral smoothing
# ==================================================
def smooth121(array_in):
    """
    Apply a 1-2-1 smoothing filter to a one-dimensional array.

    Interior weights: 1/4, 1/2, 1/4.
    Boundary weights: 3/4, 1/4 and 1/4, 3/4.

    The filter preserves the total sum for finite input data.
    Missing values propagate through interior points.
    """

    x = np.asarray(array_in, dtype=np.float64)
    y = x.copy()

    if x.ndim != 1:
        raise ValueError("Input must be a one-dimensional array.")

    if x.size < 2:
        return y

    # Interior points
    if x.size >= 3:
        y[1:-1] = (x[:-2] + 2*x[1:-1] + x[2:]) / 4.0

    # First point
    if np.isfinite(x[0]) and np.isfinite(x[1]):
        y[0] = (3*x[0] + x[1]) / 4.0

    # Last point
    if np.isfinite(x[-1]) and np.isfinite(x[-2]):
        y[-1] = (x[-2] + 3*x[-1]) / 4.0

    return y

# ==================================================
# 6. Symmetric and antisymmetric power
# ==================================================
def separate_power(power, wavefft, freqfft, maxwav4smth=27, verbose=True):
    """
    Calculate symmetric and antisymmetric power spectra.

    Power is averaged over temporal windows and latitude pairs.
    One pass of 1-2-1 smoothing is applied along positive
    frequencies for selected zonal wavenumbers.

    Returns
    -------
    psumanti, psumsym : numpy.ndarray
        Antisymmetric and symmetric power spectra.
    """

    nWindow, nfreq, nlat, nwave = power.shape
    half = nlat // 2

    # Northern Hemisphere: antisymmetric component
    psumanti = np.mean(power[:, :, :half, :], axis=(0, 2))

    # Southern Hemisphere: symmetric component
    psumsym = np.mean(power[:, :, -half:, :], axis=(0, 2))

    # Select positive frequencies and smoothing wavenumbers
    indfreq = np.flatnonzero(freqfft > 0)
    indwave = np.flatnonzero(np.abs(wavefft) <= maxwav4smth)

    # Apply one 1-2-1 smoothing pass along frequency
    for wv in indwave:
        psumanti[indfreq, wv] = smooth121(psumanti[indfreq, wv])
        psumsym[indfreq, wv] = smooth121(psumsym[indfreq, wv])

    # Remove zero-frequency power
    psumanti[freqfft == 0, :] = np.nan
    psumsym[freqfft == 0, :] = np.nan

    if verbose:
        print("\nPyWKSpectra - Symmetric and Antisymmetric Power")
        print("=" * 55)
        print(f"Spectral windows:  {nWindow}")
        print(f"Latitude pairs:    {half}")
        print("Frequency smoothing: 1 pass (1-2-1)")
        print(f"Smoothing range:   +/-{maxwav4smth} wavenumbers")
        print(f"Output shape:      {psumsym.shape}")
        print("=" * 55)

    return psumanti, psumsym

# ==================================================
# 7. Background power spectrum
# ==================================================
def derive_background(power, wavefft, freqfft, maxwav4smth=27, progress=True):
    """
    Derive the Wheeler-Kiladis background power spectrum.

    The background spectrum is obtained by averaging spectral
    power over all non-equatorial latitudes and temporal windows.

    Wavenumber smoothing:
        f < 0.1          :  5 passes
        0.1 <= f < 0.2   : 10 passes
        0.2 <= f < 0.3   : 20 passes
        f >= 0.3         : 40 passes

    Frequency smoothing:
        10 passes for 0 < f <= 0.8 cycles/day.

    Parameters
    ----------
    power : numpy.ndarray
        Spectral power (window, frequency, latitude, wavenumber).
    wavefft : numpy.ndarray
        Zonal wavenumbers.
    freqfft : numpy.ndarray
        Frequencies in cycles per day.
    maxwav4smth : int, default 27
        Maximum absolute wavenumber for smoothing.
    progress : bool, default True
        Display processing information and progress bars.

    Returns
    -------
    psumb : numpy.ndarray
        Background power spectrum (frequency, wavenumber).
    """

    nWindow, nfreq, nlat, nwave = power.shape
    half = nlat // 2

    # Average over all non-equatorial latitudes and windows
    north = np.sum(power[:, :, :half, :], axis=(0, 2), dtype=np.float64)
    south = np.sum(power[:, :, -half:, :], axis=(0, 2), dtype=np.float64)
    psumb = (north + south) / (nWindow * 2 * half)

    # Select smoothing coordinates
    indwave = np.flatnonzero(np.abs(wavefft) <= maxwav4smth)
    indfreq = np.flatnonzero(freqfft > 0)
    indfreqsmth = np.flatnonzero((freqfft > 0) & (freqfft <= 0.8))

    if progress:
        print("\nPyWKSpectra - Background Power Spectrum")
        print("=" * 55)
        print(f"Spectral windows:  {nWindow}")
        print(f"Latitude points:   {2 * half}")
        print(f"Smoothing range:   +/-{maxwav4smth} wavenumbers")
        print("-" * 55)

    # --------------------------------------------------
    # 1. Wavenumber smoothing
    # --------------------------------------------------

    for j, tt in enumerate(indfreq):

        f = freqfft[tt]

        if f < 0.1:
            passes = 5
        elif f < 0.2:
            passes = 10
        elif f < 0.3:
            passes = 20
        else:
            passes = 40

        spectrum = psumb[tt, indwave].copy()

        for i in range(passes):
            spectrum = smooth121(spectrum)

        psumb[tt, indwave] = spectrum

        if progress:
            percent = 100 * (j + 1) / len(indfreq)
            print(f"\rWavenumber smoothing: {percent:5.1f}%", end="", flush=True)

    if progress:
        print()

    # --------------------------------------------------
    # 2. Frequency smoothing
    # --------------------------------------------------

    for j, wv in enumerate(indwave):

        spectrum = psumb[indfreqsmth, wv].copy()

        for i in range(10):
            spectrum = smooth121(spectrum)

        psumb[indfreqsmth, wv] = spectrum

        if progress:
            percent = 100 * (j + 1) / len(indwave)
            print(f"\rFrequency smoothing:  {percent:5.1f}%", end="", flush=True)

    # Remove zero-frequency power
    psumb[freqfft == 0, :] = np.nan

    if progress:
        print("\nProcessing completed.")
        print(f"Output shape:       {psumb.shape}")
        print("=" * 55)

    return psumb

# ==================================================
# 8. Main Wheeler-Kiladis analysis
# ==================================================
def compute_wk_spectra(source, varname, latBound, nDayWin=96, nDaySkip=-48, max_wn=20, max_freq=0.5, verbose=True):
    """
    Compute Wheeler-Kiladis space-time power spectra.

    This function performs the complete spectral analysis, including
    data preparation, detrending, low-frequency filtering, equatorial
    decomposition, space-time Fourier transforms, spectral averaging,
    background estimation, and power normalization.

    Parameters
    ----------
    source : str, Path, xr.Dataset, or xr.DataArray
        Input atmospheric data (NetCDF, GRIB, or xarray).
    varname : str
        Name of the atmospheric variable to analyze.
    latBound : float
        Maximum absolute latitude for the tropical domain.
        Automatically adjusted to available symmetric coverage.
    nDayWin : int, default 96
        Length of each spectral window in days.
    nDaySkip : int, default -48
        Days skipped between consecutive windows.
        Negative values indicate overlapping windows.
    max_wn : int, default 20
        Maximum absolute zonal wavenumber in the output.
    max_freq : float, default 0.5
        Maximum frequency in cycles per day.
    verbose : bool, default True
        Print processing information and progress bars.

    Returns
    -------
    results : xarray.Dataset
        Dataset containing symmetric, antisymmetric, background,
        logarithmic, and background-normalized power spectra.

        Dimensions:
            frequency, wavenumber

    Notes
    -----
    The analysis follows the Wheeler-Kiladis methodology, using
    a Hann temporal taper and the spectral processing conventions
    implemented in PyWKSpectra.
    """

    # --------------------------------------------------
    # 1. Load and validate input data
    # --------------------------------------------------

    if not isinstance(varname, str) or not varname.strip():
        raise ValueError("varname must be a non-empty string.")

    if isinstance(source, xr.DataArray) and source.name not in (None, varname):
        raise ValueError(f"DataArray name {source.name!r} does not match varname {varname!r}.")

    data, spd = load_wk_data(source, varname=varname, latBound=latBound, verbose=verbose)

    if data.name is None:
        data = data.rename(varname)

    # --------------------------------------------------
    # 2. Calculate temporal sampling parameters
    # --------------------------------------------------

    ntim, nlat, nlon, nDayTot, nSampWin, nSampSkip = sampling_vars(data, spd, nDayWin, nDaySkip, verbose=verbose)

    # --------------------------------------------------
    # 3. Remove dominant low-frequency signals
    # --------------------------------------------------

    array_dt = remove_dominant_signals(data, spd, nDayTot, nSampWin, progress=verbose)

    # --------------------------------------------------
    # 4. Equatorial symmetric/antisymmetric decomposition
    # --------------------------------------------------

    array_as = decompose_symasym(array_dt, verbose=verbose)

    # Release intermediate data
    del array_dt

    # --------------------------------------------------
    # 5. Compute space-time Fourier coefficients
    # --------------------------------------------------

    wavefft, freqfft, peeAS = spectral_coefficients(array_as, spd, nSampWin, nSampSkip, progress=verbose)

    del array_as

    # --------------------------------------------------
    # 6. Calculate spectral power
    # --------------------------------------------------

    power = np.abs(peeAS)**2

    del peeAS

    # --------------------------------------------------
    # 7. Separate symmetric and antisymmetric power
    # --------------------------------------------------

    psumanti, psumsym = separate_power(power, wavefft, freqfft, verbose=verbose)

    # --------------------------------------------------
    # 8. Derive background spectrum
    # --------------------------------------------------

    psumb = derive_background(power, wavefft, freqfft, progress=verbose)

    del power

    # --------------------------------------------------
    # 9. Crop spectral domain
    # --------------------------------------------------

    indwave = np.flatnonzero(np.abs(wavefft) <= max_wn)
    indfreq = np.flatnonzero((freqfft > 0) & (freqfft <= max_freq))

    if len(indwave) == 0 or len(indfreq) == 0:
        raise ValueError("The requested spectral domain contains no available frequencies or wavenumbers.")

    wavefft = wavefft[indwave]
    freqfft = freqfft[indfreq]

    psumanti = psumanti[np.ix_(indfreq, indwave)]
    psumsym = psumsym[np.ix_(indfreq, indwave)]
    psumb = psumb[np.ix_(indfreq, indwave)]

    # --------------------------------------------------
    # 10. Validate spectral power
    # --------------------------------------------------

    if not np.isfinite(psumanti).all() or not np.isfinite(psumsym).all():
        raise ValueError("Invalid symmetric or antisymmetric spectral power values detected.")

    if not np.isfinite(psumb).all() or np.any(psumb <= 0):
        raise ValueError("Background power must contain strictly positive, finite values.")

    if np.any(psumanti < 0) or np.any(psumsym < 0):
        raise ValueError("Spectral power must be non-negative.")

    # --------------------------------------------------
    # 11. Calculate logarithmic spectra
    # --------------------------------------------------

    with np.errstate(divide="ignore"):
        psumanti_log = np.log10(psumanti)
        psumsym_log = np.log10(psumsym)
        psumb_log = np.log10(psumb)

    psumanti_log[~np.isfinite(psumanti_log)] = np.nan
    psumsym_log[~np.isfinite(psumsym_log)] = np.nan

    # --------------------------------------------------
    # 12. Normalize spectra by background power
    # --------------------------------------------------

    psumanti_r = psumanti / psumb
    psumsym_r = psumsym / psumb

    # --------------------------------------------------
    # 13. Create xarray Dataset
    # --------------------------------------------------

    coords = {
        "frequency": ("frequency", freqfft, {"long_name": "Frequency", "units": "cycles/day"}),
        "wavenumber": ("wavenumber", wavefft, {"long_name": "Zonal wavenumber", "units": "1"})
    }

    results = xr.Dataset(
        data_vars={
            "antisymmetric_power": (("frequency", "wavenumber"), psumanti),
            "symmetric_power": (("frequency", "wavenumber"), psumsym),
            "background_power": (("frequency", "wavenumber"), psumb),
            "antisymmetric_log_power": (("frequency", "wavenumber"), psumanti_log),
            "symmetric_log_power": (("frequency", "wavenumber"), psumsym_log),
            "background_log_power": (("frequency", "wavenumber"), psumb_log),
            "antisymmetric_normalized": (("frequency", "wavenumber"), psumanti_r),
            "symmetric_normalized": (("frequency", "wavenumber"), psumsym_r)
        },
        coords=coords
    )

    # --------------------------------------------------
    # 14. Add global metadata
    # --------------------------------------------------

    nWindow = 1 + (ntim - nSampWin) // (nSampWin + nSampSkip)

    results.attrs = {
        "title": "Wheeler-Kiladis Space-Time Spectra",
        "software": "PyWKSpectra",
        "author": "Alejandro Jaramillo Moreno",
        "institution": "Instituto de Ciencias de la Atmosfera y Cambio Climatico, UNAM",
        "variable": data.name,
        "start_date": str(data.time.values[0])[:10],
        "end_date": str(data.time.values[-1])[:10],
        "latitude_bound": float(np.max(np.abs(data.lat.values))),
        "samples_per_day": int(spd),
        "window_days": int(nDayWin),
        "window_skip_days": int(nDaySkip),
        "number_of_windows": int(nWindow),
        "temporal_taper": "Hann",
        "low_frequency_cutoff": float(1.0 / nDayWin),
        "equator_in_latitude_average": "Excluded",
        "method": "Wheeler-Kiladis space-time spectral analysis"
    }

    # --------------------------------------------------
    # 15. Add variable attributes
    # --------------------------------------------------

    results["antisymmetric_power"].attrs = {"long_name": "Antisymmetric spectral power"}
    results["symmetric_power"].attrs = {"long_name": "Symmetric spectral power"}
    results["background_power"].attrs = {"long_name": "Background spectral power"}

    results["antisymmetric_log_power"].attrs = {"long_name": "Log10 antisymmetric spectral power"}
    results["symmetric_log_power"].attrs = {"long_name": "Log10 symmetric spectral power"}
    results["background_log_power"].attrs = {"long_name": "Log10 background spectral power"}

    results["antisymmetric_normalized"].attrs = {"long_name": "Background-normalized antisymmetric power", "units": "1"}
    results["symmetric_normalized"].attrs = {"long_name": "Background-normalized symmetric power", "units": "1"}

    # --------------------------------------------------
    # 16. Final summary
    # --------------------------------------------------

    if verbose:
        print("\nPyWKSpectra - Analysis Completed")
        print("=" * 55)
        print(f"Variable:           {data.name}")
        print(f"Analysis period:    {results.attrs['start_date']} to {results.attrs['end_date']}")
        print(f"Latitude bound:     +/-{results.attrs['latitude_bound']:g} degrees")
        print(f"Spectral windows:   {nWindow}")
        print(f"Window length:      {nDayWin} days")
        print(f"Window overlap:     {max(0, -nDaySkip)} days")
        print(f"Wavenumber range:   {wavefft.min():g} to {wavefft.max():g}")
        print(f"Frequency range:    {freqfft.min():.4f} to {freqfft.max():.4f} cycles/day")
        print(f"Output dimensions:  {psumsym.shape}")
        print("Output format:      xarray.Dataset")
        print("=" * 55)

    return results