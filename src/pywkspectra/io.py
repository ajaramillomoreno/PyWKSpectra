# Copyright (C) 2026 Alejandro Jaramillo-Moreno
# SPDX-License-Identifier: GPL-3.0-only
"""
PyWKSpectra - Data Input and Validation
=======================================

Wheeler-Kiladis Space-Time Spectral Analysis in Python

Author: Alejandro Jaramillo-Moreno
        Departamento de Cambio Climático
        Instituto de Ciencias de la Atmósfera y Cambio Climático
        Universidad Nacional Autónoma de México (UNAM)

Description
-----------
Functions for loading, organizing, and validating atmospheric
data before Wheeler-Kiladis spectral analysis.
"""

from pathlib import Path
import numpy as np
import xarray as xr


def load_wk_data(source, varname=None, latBound=15, verbose=True):
    """
    Load atmospheric data and prepare them for Wheeler-Kiladis analysis.

    Parameters
    ----------
    source : str, Path, xr.Dataset, or xr.DataArray
        Input NetCDF/GRIB file or xarray object.
    varname : str, optional
        Variable name when the input is a Dataset.
    latBound : float, default 15
        Maximum absolute latitude for the tropical domain.
        Automatically reduced if symmetric coverage is unavailable.
    verbose : bool, default True
        Print a summary of the imported data.

    Returns
    -------
    data : xr.DataArray
        Validated data with dimensions (time, lat, lon).
        Latitudes are descending and symmetric about the equator.
        Longitudes are ascending within [-180, 180).
    spd : int
        Number of samples per day.

    Notes
    -----
    Supports rectilinear grids and datetime64 time coordinates.
    GRIB files require the optional cfgrib/ecCodes backend.
    """

    # --------------------------------------------------
    # 1. Validate parameters
    # --------------------------------------------------

    if isinstance(latBound, (bool, np.bool_)) or not isinstance(latBound, (int, float, np.number)) or not np.isfinite(latBound) or not 0 < latBound <= 90:
        raise ValueError("latBound must be a finite number between 0 and 90 degrees.")

    opened = None

    # --------------------------------------------------
    # 2. Open input data
    # --------------------------------------------------

    if isinstance(source, xr.DataArray):
        data = source.copy(deep=False)

        if varname is not None and data.name not in (None, varname):
            raise ValueError(f"DataArray name {data.name!r} does not match varname {varname!r}.")

    else:
        if isinstance(source, xr.Dataset):
            ds = source

        else:
            path = Path(source)

            if not path.is_file():
                raise FileNotFoundError(f"File not found: {path}")

            is_grib = path.suffix.lower() in (".grib", ".grib2", ".grb", ".grb2")

            try:
                ds = xr.open_dataset(path, engine="cfgrib" if is_grib else None)

            except ValueError as exc:
                if is_grib:
                    raise ImportError("Cannot open GRIB input. Install cfgrib and ecCodes, and verify the GRIB structure.") from exc
                raise

            opened = ds

    try:

        # --------------------------------------------------
        # 3. Select variable
        # --------------------------------------------------

        if not isinstance(source, xr.DataArray):

            if varname is None:
                if len(ds.data_vars) != 1:
                    raise ValueError(f"Specify varname. Available variables: {list(ds.data_vars)}")

                varname = next(iter(ds.data_vars))

            if varname not in ds.data_vars:
                raise ValueError(f"Variable {varname!r} not found. Available variables: {list(ds.data_vars)}")

            data = ds[varname]

        elif data.name is None and varname is not None:
            data = data.rename(varname)

        # --------------------------------------------------
        # 4. Identify and standardize coordinates
        # --------------------------------------------------

        options = {"time": ("time", "valid_time", "date"), "lat": ("lat", "latitude"), "lon": ("lon", "longitude")}
        standards = {"time": "time", "lat": "latitude", "lon": "longitude"}

        for standard, candidates in options.items():

            matching = [name for name in candidates if name in data.coords and data[name].ndim == 1]

            if not matching:
                matching = [name for name in data.coords if data[name].ndim == 1 and str(data[name].attrs.get("standard_name", "")).lower() == standards[standard]]

            if len(matching) != 1:
                raise ValueError(f"Cannot identify a unique {standard} coordinate. Found: {list(data.coords)}")

            old_name = matching[0]
            old_dim = data[old_name].dims[0]

            if old_dim != old_name:
                data = data.swap_dims({old_dim: old_name})

            if old_name != standard:
                data = data.rename({old_name: standard})

        if set(data.dims) != {"time", "lat", "lon"} or data.ndim != 3:
            raise ValueError(f"Expected dimensions (time, lat, lon); found {data.dims}.")

        data = data.transpose("time", "lat", "lon")

        # --------------------------------------------------
        # 5. Validate time and calculate SPD
        # --------------------------------------------------

        time = np.asarray(data.time.values)

        if not np.issubdtype(time.dtype, np.datetime64) or len(time) < 2 or np.isnat(time).any():
            raise ValueError("Time must contain at least two valid datetime64 values.")

        steps = np.diff(time.astype("datetime64[ns]")).astype("timedelta64[ns]").astype(np.int64)

        if steps[0] <= 0 or not np.all(steps == steps[0]):
            raise ValueError("Time must be strictly increasing, evenly spaced, and have no missing timestamps.")

        ns_per_day = 86_400_000_000_000

        if ns_per_day % int(steps[0]) != 0:
            raise ValueError("The sampling interval must divide one day exactly.")

        spd = int(ns_per_day // steps[0])

        # --------------------------------------------------
        # 6. Select symmetric tropical latitudes
        # --------------------------------------------------

        latitude = np.asarray(data.lat.values, dtype=float)

        if not np.isfinite(latitude).all() or not np.any(latitude > 0) or not np.any(latitude < 0):
            raise ValueError("Input must contain valid latitudes in both hemispheres.")

        available_bound = min(float(latBound), float(latitude[latitude > 0].max()), float(-latitude[latitude < 0].min()))

        data = data.sortby("lat", ascending=False)
        data = data.sel(lat=slice(available_bound, -available_bound))

        selected_lat = np.asarray(data.lat.values, dtype=float)

        if selected_lat.size < 2 or not np.allclose(selected_lat, -selected_lat[::-1], rtol=0, atol=1e-5):
            raise ValueError("No valid symmetric latitude grid within latBound; opposite latitude pairs are required.")

        if selected_lat.size > 2 and not np.allclose(np.diff(selected_lat), np.diff(selected_lat)[0], rtol=1e-5, atol=1e-5):
            raise ValueError("Selected latitudes must be regularly spaced.")

        actual_bound = float(np.max(np.abs(selected_lat)))

        # --------------------------------------------------
        # 7. Normalize and validate longitudes
        # --------------------------------------------------

        longitude = np.asarray(data.lon.values, dtype=float)

        if not np.isfinite(longitude).all():
            raise ValueError("Longitude coordinates must be finite.")

        data = data.assign_coords(lon=(longitude + 180) % 360 - 180)
        data = data.sortby("lon")

        lon = np.asarray(data.lon.values, dtype=float)

        if lon.size < 2 or not np.all(np.diff(lon) > 0):
            raise ValueError("Longitude coordinates must be unique.")

        lon_step = np.diff(lon)

        if not np.allclose(lon_step, lon_step[0], rtol=1e-5, atol=1e-5) or not np.isclose(lon.size * lon_step[0], 360, rtol=1e-5, atol=1e-4):
            raise ValueError("Wheeler-Kiladis analysis requires a regular, complete 360-degree longitude grid.")

        # --------------------------------------------------
        # 8. Validate missing values
        # --------------------------------------------------

        data = data.load()

        if not np.issubdtype(data.dtype, np.number) or np.issubdtype(data.dtype, np.complexfloating) or not np.isfinite(data.values).all():
            raise ValueError("Input field must be real-valued and contain no missing or non-finite values.")

        # --------------------------------------------------
        # 9. Print data summary
        # --------------------------------------------------

        if verbose:

            print("\nPyWKSpectra - Input Data")
            print("=" * 55)

            print(f"Variable:           {data.name or 'unnamed'}")
            print(f"Start date:         {str(data.time.values[0])[:19]}")
            print(f"End date:           {str(data.time.values[-1])[:19]}")
            print(f"Requested latBound: +/-{latBound:g} degrees")
            print(f"Applied latBound:   +/-{actual_bound:g} degrees")

            if actual_bound < float(latBound) - 1e-5:
                print(f"NOTE: latBound adjusted from +/-{latBound:g} to +/-{actual_bound:g} degrees for symmetric coverage.")

            print(f"Latitude range:     {selected_lat[-1]:g} to {selected_lat[0]:g} degrees")
            print(f"Longitude range:    {lon[0]:g} to {lon[-1]:g} degrees")
            print(f"Samples per day:    {spd}")
            print(f"Data shape:         {data.shape}")
            print("Missing values:     None")

            print("=" * 55)

        return data, spd

    finally:
        if opened is not None:
            opened.close()
