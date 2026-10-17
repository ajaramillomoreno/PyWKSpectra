
"""
PyWKSpectra - Spectral Visualization
====================================

Wheeler-Kiladis Space-Time Spectral Analysis in Python

Author: Alejandro Jaramillo-Moreno
        Departamento de Cambio Climático
        Instituto de Ciencias de la Atmósfera y Cambio Climático
        Universidad Nacional Autónoma de México (UNAM)

Description
-----------
Visualization functions for Wheeler-Kiladis space-time spectra
and theoretical Matsuno equatorial-wave dispersion curves.

Includes raw power spectra, background-normalized spectra,
and theoretical dispersion diagrams.
"""

from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

# ==================================================
# 1. PyWKSpectra figure signature
# ==================================================
def _add_pywkspectra_mark(ax, text_size=11):
    """Add a small software credit inside the bottom-right corner of an axis."""
    ax.text(0.985, 0.025, "PyWKSpectra", transform=ax.transAxes, ha="right", va="bottom", fontsize=max(8, text_size-3), color="#243449", fontweight="bold", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5}, zorder=30)

# ==================================================
# 2. Matsuno dispersion relations
# ==================================================
def matsuno_dispersion(he, wavenumber, meridional_modes=(1,)):
    """
    Calculate equatorial-wave dispersion curves following Matsuno (1966).

    Parameters
    ----------
    he : float
        Equivalent depth in meters.
    wavenumber : array-like
        Global zonal wavenumbers.
    meridional_modes : sequence of int
        Meridional mode numbers (n >= 1).

    Returns
    -------
    curves : dict
        Dispersion frequencies in cycles per day.
    """

    # Physical constants
    radius = 6.371008e6
    gravity = 9.80665
    rotation = 7.292e-5

    beta = 2 * rotation / radius
    c = np.sqrt(gravity * he)
    wn = np.asarray(wavenumber, dtype=float)
    k = wn / radius
    to_cpd = 86400 / (2 * np.pi)

    curves = {}

    # Kelvin wave
    curves["Kelvin"] = np.where(wn >= 0, c * k * to_cpd, np.nan)

    # Continuous MRG / EIG(n=0) branch
    curves["MRG_EIG0"] = 0.5 * (c*k + np.sqrt((c*k)**2 + 4*beta*c)) * to_cpd

    # Rossby and inertia-gravity modes (n >= 1)
    for n in meridional_modes:

        A = (c*k)**2 + (2*n + 1)*beta*c
        B = beta*c*c*k

        theta = np.arccos(np.clip(3*np.sqrt(3)*B / (2*A**1.5), -1, 1)) / 3

        gravity_branch = 2*np.sqrt(A/3) * np.cos(theta) * to_cpd
        rossby_branch = 2*np.sqrt(A/3) * np.cos(theta - 2*np.pi/3) * to_cpd

        curves[f"EIG_{n}"] = np.where(wn >= 0, gravity_branch, np.nan)
        curves[f"WIG_{n}"] = np.where(wn <= 0, gravity_branch, np.nan)
        curves[f"ER_{n}"] = np.where(wn <= 0, rossby_branch, np.nan)

    return curves

# ==================================================
# 3. Raw power spectra
# ==================================================
def plot_raw_spectra(results, figsize=(8, 5.5), text_size=11, contour_range=None, cmap="magma", max_freq_plot=None, max_wn_plot=None, freq_lines=True, cpd_lines=(3, 6, 30), dpi=120, save_dir=None, show=True):
    """
    Plot raw Wheeler-Kiladis power spectra.

    Generates three publication-style figures:
    antisymmetric, symmetric, and background log10 power spectra.

    Parameters
    ----------
    results : xarray.Dataset
        Output from compute_wk_spectra().
    figsize : tuple, default (8, 5.5)
        Figure dimensions in inches.
    text_size : int, default 11
        Base font size.
    contour_range : tuple or None, default None
        Color levels (minimum, maximum, interval).
        If None, shared limits are calculated automatically.
    cmap : str, default 'magma'
        Matplotlib colormap.
    max_freq_plot : float, optional
        Maximum plotted frequency in cycles/day.
    max_wn_plot : float, optional
        Maximum absolute zonal wavenumber.
    freq_lines : bool, default True
        Display reference period lines.
    cpd_lines : sequence, default (3, 6, 30)
        Periods in days for reference lines.
    dpi : int, default 120
        Resolution for figures and PNG export.
    save_dir : str or Path, optional
        Directory for saving PNG and PDF figures.
    show : bool, default True
        Display figures.

    Returns
    -------
    figures : dict
        Dictionary of Matplotlib Figure objects.
    """

    # --------------------------------------------------
    # 1. Extract spectra and coordinates
    # --------------------------------------------------

    fields = {
        "antisymmetric_raw_log": "antisymmetric_log_power",
        "symmetric_raw_log": "symmetric_log_power",
        "background_log": "background_log_power"
    }

    if not isinstance(results, xr.Dataset):
        raise TypeError("results must be an xarray.Dataset.")

    missing = [name for name in fields.values() if name not in results]

    if missing:
        raise ValueError(f"Missing spectral variables: {missing}")

    wave = np.asarray(results.wavenumber.values, dtype=float)
    freq = np.asarray(results.frequency.values, dtype=float)

    # --------------------------------------------------
    # 2. Select plotting domain
    # --------------------------------------------------

    available_wn = float(np.max(np.abs(wave)))
    available_freq = float(np.max(freq))

    max_wn_plot = available_wn if max_wn_plot is None else float(max_wn_plot)
    max_freq_plot = available_freq if max_freq_plot is None else float(max_freq_plot)

    if max_wn_plot <= 0 or max_freq_plot <= 0:
        raise ValueError("Plot limits must be positive.")

    if max_wn_plot > available_wn:
        print(f"Requested wavenumber exceeds available data. Using +/-{available_wn:g}.")
        max_wn_plot = available_wn

    if max_freq_plot > available_freq:
        print(f"Requested frequency exceeds available data. Using {available_freq:g} cycles/day.")
        max_freq_plot = available_freq

    indwave = np.flatnonzero(np.abs(wave) <= max_wn_plot)
    indfreq = np.flatnonzero((freq > 0) & (freq <= max_freq_plot))

    if len(indwave) < 2 or len(indfreq) < 2:
        raise ValueError("Selected domain must include at least two frequencies and wavenumbers.")

    waveplt = wave[indwave]
    freqplt = freq[indfreq]

    spectra = {key: results[name].isel(frequency=indfreq, wavenumber=indwave).to_numpy() for key, name in fields.items()}

    # --------------------------------------------------
    # 3. Define shared color levels
    # --------------------------------------------------

    finite = np.concatenate([values[np.isfinite(values)] for values in spectra.values()])

    if finite.size == 0:
        raise ValueError("No finite spectral power values in the selected domain.")

    if contour_range is None:
        vmin, vmax = np.percentile(finite, [2, 98])
        padding = max(0.05 * (vmax - vmin), 0.01)
        levels = np.linspace(vmin - padding, vmax + padding, 19)

    else:
        vmin, vmax, step = map(float, contour_range)

        if not (np.isfinite([vmin, vmax, step]).all() and vmin < vmax and step > 0):
            raise ValueError("contour_range must contain valid minimum, maximum, and positive interval.")

        levels = np.arange(vmin, vmax + step * 0.5, step)

        if len(levels) < 2:
            raise ValueError("At least two contour levels are required.")

    # --------------------------------------------------
    # 4. Figure titles and metadata
    # --------------------------------------------------

    titles = {
        "antisymmetric_raw_log": "Antisymmetric power spectrum",
        "symmetric_raw_log": "Symmetric power spectrum",
        "background_log": "Background power spectrum"
    }

    metadata = []

    if results.attrs.get("variable"):
        metadata.append(str(results.attrs["variable"]).upper())

    if results.attrs.get("start_date") and results.attrs.get("end_date"):
        metadata.append(f"{results.attrs['start_date']} to {results.attrs['end_date']}")

    if results.attrs.get("latitude_bound") is not None:
        lat = float(results.attrs["latitude_bound"])
        metadata.append(f"{lat:g}°S–{lat:g}°N")

    subtitle = "  |  ".join(metadata)

    # --------------------------------------------------
    # 5. Publication-style settings
    # --------------------------------------------------

    style = {
        "font.family": "DejaVu Sans",
        "font.size": text_size,
        "axes.labelsize": text_size,
        "axes.titlesize": text_size + 1,
        "xtick.labelsize": text_size - 1,
        "ytick.labelsize": text_size - 1,
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white"
    }

    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

    figures = {}

    # --------------------------------------------------
    # 6. Generate figures
    # --------------------------------------------------

    with plt.rc_context(style):

        for key, values in spectra.items():

            fig, ax = plt.subplots(figsize=figsize, dpi=dpi, layout="constrained")

            # Filled contours
            cf = ax.contourf(waveplt, freqplt, np.ma.masked_invalid(values), levels=levels, cmap=cmap, extend="both")

            # Subtle contour lines
            ax.contour(waveplt, freqplt, np.ma.masked_invalid(values), levels=levels[::4], colors="#1a2230", linewidths=0.3, alpha=0.30)

            # Reference period lines
            if freq_lines:

                for period in cpd_lines:

                    if period <= 0:
                        raise ValueError("Reference periods must be positive.")

                    fline = 1.0 / period

                    if freqplt[0] <= fline <= max_freq_plot:

                        ax.axhline(fline, color="white", linewidth=0.8, linestyle=(0, (3, 3)), alpha=0.85)

                        ax.text(0.985, fline + 0.005 * max_freq_plot, f"{period:g} d",
                                transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                                fontsize=text_size-2, color="white",
                                bbox={"facecolor": "#1f2937", "edgecolor": "none", "alpha": 0.6, "pad": 1.3})

            # Zero wavenumber
            ax.axvline(0, color="white", linewidth=0.65, alpha=0.70)

            # Axes
            ax.set_xlim(-max_wn_plot, max_wn_plot)
            ax.set_ylim(float(freqplt[0]), max_freq_plot)
            ax.set_xlabel("Zonal wavenumber")
            ax.set_ylabel("Frequency (cycles day$^{-1}$)")

            # Title and subtitle
            ax.set_title(titles[key], loc="left", fontweight="bold", pad=34)

            if subtitle:
                ax.text(0, 1.02, subtitle, transform=ax.transAxes,
                        ha="left", va="bottom", fontsize=text_size-2, color="#667085")

            # Axis styling
            ax.tick_params(direction="out", length=3.5, width=0.8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # Colorbar
            cbar = fig.colorbar(cf, ax=ax, fraction=0.047, pad=0.035, shrink=0.93, format="%.2f")
            cbar.set_label(r"$\log_{10}$(spectral power)", labelpad=10)
            cbar.set_ticks(np.linspace(levels[0], levels[-1], 6))
            cbar.ax.tick_params(length=3, width=0.7)
            cbar.outline.set_linewidth(0.7)

            # Discreet software credit inside the spectrum
            _add_pywkspectra_mark(ax, text_size=text_size)

            figures[key] = fig

            # Save figures
            if save_dir is not None:
                fig.savefig(save_dir / f"wk_{key}.png", dpi=dpi, bbox_inches="tight")
                fig.savefig(save_dir / f"wk_{key}.pdf", bbox_inches="tight")

        if show:
            plt.show()

    return figures

# ==================================================
# 4. Background-normalized WK spectra
# ==================================================
def plot_wk_spectra(results, figsize=(12, 5.4), dpi=120, text_size=11, contour_range=(0, 2, 0.1), cmap="RdBu_r", max_wn_plot=None, max_freq_plot=None, matsuno_lines=True, he=(12, 25, 50), meridional_modes=(1,), freq_lines=True, cpd_lines=(3, 6, 30), contour_lines=True, labels=False, save_dir=None, show=True):
    """
    Plot symmetric and antisymmetric Wheeler-Kiladis spectra.

    Generates a two-panel figure using background-normalized
    power spectra and optional theoretical Matsuno dispersion curves.

    Parameters
    ----------
    results : xarray.Dataset
        Output from compute_wk_spectra().
    figsize : tuple, default (12, 5.4)
        Figure dimensions in inches.
    dpi : int, default 120
        Figure resolution.
    text_size : int, default 11
        Base font size.
    contour_range : tuple, default (0, 2, 0.1)
        Contour minimum, maximum, and interval.
    cmap : str, default 'RdBu_r'
        Matplotlib colormap.
    max_wn_plot : float, optional
        Maximum absolute plotted wavenumber.
    max_freq_plot : float, optional
        Maximum plotted frequency in cycles/day.
    matsuno_lines : bool, default True
        Overlay theoretical dispersion curves.
    he : sequence, default (12, 25, 50)
        Equivalent depths in meters.
    meridional_modes : sequence, default (1,)
        Meridional modes for Rossby and inertia-gravity waves.
    freq_lines : bool, default True
        Draw reference-period lines.
    cpd_lines : sequence, default (3, 6, 30)
        Reference periods in days.
    contour_lines : bool, default True
        Draw selected contour lines above normalized power 1.
    labels : bool, default False
        Label selected wave branches.
    save_dir : str or Path, optional
        Directory for PNG and PDF output.
    show : bool, default True
        Display figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray
        Axes for symmetric and antisymmetric panels.
    """

    # --------------------------------------------------
    # 1. Extract data and select plotting domain
    # --------------------------------------------------

    fields = ("symmetric_normalized", "antisymmetric_normalized")

    if not isinstance(results, xr.Dataset) or any(name not in results for name in fields):
        raise ValueError(f"results must contain {fields}.")

    wave = np.asarray(results.wavenumber.values, dtype=float)
    freq = np.asarray(results.frequency.values, dtype=float)

    available_wn = np.max(np.abs(wave))
    available_freq = np.max(freq)

    max_wn_plot = available_wn if max_wn_plot is None else min(float(max_wn_plot), available_wn)
    max_freq_plot = available_freq if max_freq_plot is None else min(float(max_freq_plot), available_freq)

    if max_wn_plot <= 0 or max_freq_plot <= 0:
        raise ValueError("Plot limits must be positive.")

    sub = results.sel(wavenumber=slice(-max_wn_plot, max_wn_plot), frequency=slice(0, max_freq_plot))

    w = sub.wavenumber.values
    f = sub.frequency.values

    if len(w) < 2 or len(f) < 2:
        raise ValueError("Selected spectral domain contains insufficient points.")

    vmin, vmax, interval = contour_range

    if vmin >= vmax or interval <= 0:
        raise ValueError("Invalid contour range.")

    levels = np.arange(vmin, vmax + interval*0.5, interval)

    # --------------------------------------------------
    # 2. Matsuno configuration
    # --------------------------------------------------

    colors = ("#243449", "#566477", "#8997a9")
    styles = ((0, (3, 2)), (0, (5, 3)), (0, (1, 2)))

    if matsuno_lines:

        if not he or any(not np.isfinite(h) or h <= 0 for h in he):
            raise ValueError("Equivalent depths must be positive.")

        if any(int(n) != n or n < 1 for n in meridional_modes):
            raise ValueError("Meridional modes must be positive integers.")

        xdisp = np.linspace(-max_wn_plot, max_wn_plot, 501)

    # --------------------------------------------------
    # 3. Figure style
    # --------------------------------------------------

    style = {
        "font.family": "DejaVu Sans",
        "font.size": text_size,
        "axes.labelsize": text_size,
        "xtick.labelsize": text_size - 1,
        "ytick.labelsize": text_size - 1,
        "axes.linewidth": 0.8
    }

    with plt.rc_context(style):

        fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi, sharex=True, sharey=True, layout="constrained")

        # --------------------------------------------------
        # 4. Plot normalized spectra
        # --------------------------------------------------

        for ax, field, title in zip(axes, fields, ("(a) Symmetric", "(b) Antisymmetric")):

            z = sub[field].transpose("frequency", "wavenumber").values

            cf = ax.contourf(w, f, z, levels=levels, cmap=cmap, extend="both")

            if contour_lines:
                line_levels = levels[(levels > 1) & (levels < vmax)][::2]
                line_levels = line_levels[(line_levels >= np.nanmin(z)) & (line_levels <= np.nanmax(z))]

                if len(line_levels) > 0:
                    ax.contour(w, f, z, levels=line_levels, colors="#263547", linewidths=0.4, alpha=0.45)

            # Zero wavenumber
            ax.axvline(0, color="#475569", linewidth=0.7, alpha=0.75)

            # Reference periods
            if freq_lines:
                for days in cpd_lines:

                    if days <= 0:
                        raise ValueError("Reference periods must be positive.")

                    if 1/days <= max_freq_plot:

                        ax.axhline(1/days, color="#3b4551", linewidth=0.7, linestyle=(0, (2, 3)), alpha=0.5)

                        if ax is axes[0]:
                            ax.annotate(f"{days:g} d", (w.min(), 1/days), xytext=(4, 3), textcoords="offset points", fontsize=text_size-3, color="#333d4b", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1})

            ax.set_title(title, loc="left", fontweight="bold", pad=10)
            ax.set_xlim(-max_wn_plot, max_wn_plot)
            ax.set_ylim(float(f[0]), max_freq_plot)
            ax.set_xlabel("Zonal wavenumber")

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(direction="out", length=3.5)

        axes[0].set_ylabel("Frequency (cycles day$^{-1}$)")

        # --------------------------------------------------
        # 5. Overlay Matsuno dispersion curves
        # --------------------------------------------------

        if matsuno_lines:

            for j, depth in enumerate(he):

                curves = matsuno_dispersion(depth, xdisp, meridional_modes)

                line_style = {
                    "color": colors[j % len(colors)],
                    "linestyle": styles[j % len(styles)],
                    "linewidth": 1.05,
                    "alpha": 0.9,
                    "zorder": 5
                }

                # Symmetric modes
                axes[0].plot(xdisp, curves["Kelvin"], **line_style)

                # Antisymmetric n=0 branch
                axes[1].plot(xdisp, curves["MRG_EIG0"], **line_style)

                # Higher meridional modes
                for n in meridional_modes:

                    ax = axes[0] if n % 2 == 1 else axes[1]

                    for mode in (f"ER_{n}", f"EIG_{n}", f"WIG_{n}"):
                        ax.plot(xdisp, curves[mode], **line_style)

            # Equivalent depth legend
            handles = [
                Line2D([0], [0], color=colors[j % len(colors)], linestyle=styles[j % len(styles)], linewidth=1.4, label=f"{depth:g} m")
                for j, depth in enumerate(he)
            ]

            axes[1].legend(handles=handles, title="Equivalent depth", loc="upper right", frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, fontsize=text_size-3, title_fontsize=text_size-3)

            # Optional wave labels
            if labels:

                curves = matsuno_dispersion(he[len(he)//2], xdisp, meridional_modes)

                positions = [
                    (axes[0], "Kelvin", 0.60*max_wn_plot, "Kelvin"),
                    (axes[1], "MRG_EIG0", -0.40*max_wn_plot, "MRG"),
                    (axes[1], "MRG_EIG0", 0.25*max_wn_plot, "EIG n=0")
                ]

                for ax, key, x, label in positions:

                    y = float(np.interp(x, xdisp, curves[key]))

                    if 0 < y < 0.9*max_freq_plot:
                        ax.text(x, y + 0.015*max_freq_plot, label, ha="center", fontsize=text_size-2, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.2})

        # --------------------------------------------------
        # 6. Shared colorbar and title
        # --------------------------------------------------

        cbar = fig.colorbar(cf, ax=axes, location="bottom", fraction=0.06, pad=0.12, shrink=0.72, aspect=36)
        cbar.set_label("Power / background power")
        cbar.set_ticks(np.linspace(vmin, vmax, 5))

        fig.suptitle("Wheeler–Kiladis spectra", fontsize=text_size+3, fontweight="bold")

        # Credit on the right panel, away from the shared colorbar
        _add_pywkspectra_mark(axes[1], text_size=text_size)

        # --------------------------------------------------
        # 7. Save and display
        # --------------------------------------------------

        if save_dir is not None:

            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)

            fig.savefig(save_dir / "wk_normalized_spectra.png", dpi=dpi, bbox_inches="tight")
            fig.savefig(save_dir / "wk_normalized_spectra.pdf", bbox_inches="tight")

        if show:
            plt.show()

    return fig, axes


# ==================================================
# 5. Matsuno dispersion diagrams
# ==================================================
def plot_matsuno(he=25, max_wn=20, max_freq=0.5, meridional_modes=(1, 2, 3), figsize=(9, 6), dpi=120, text_size=11, save_path=None):
    """
    Plot theoretical equatorial-wave dispersion curves following Matsuno (1966).

    Parameters
    ----------
    he : float, default 25
        Equivalent depth in meters.
    max_wn : int, default 20
        Maximum absolute zonal wavenumber.
    max_freq : float, default 0.5
        Maximum frequency in cycles per day.
    meridional_modes : tuple, default (1, 2, 3)
        Meridional mode numbers.
    figsize : tuple, default (9, 6)
        Figure size in inches.
    dpi : int, default 120
        Figure resolution.
    text_size : int, default 11
        Base font size.
    save_path : str or Path, optional
        Path for saving the figure.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
        Figure and axes objects.
    """

    # Calculate Matsuno dispersion curves
    wn = np.linspace(-max_wn, max_wn, 1001)
    curves = matsuno_dispersion(he, wn, meridional_modes)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Kelvin wave
    ax.plot(wn, curves["Kelvin"], color="#D55E00", lw=2.0, label="Kelvin")

    # Mixed Rossby-Gravity / EIG(n=0)
    ax.plot(wn, curves["MRG_EIG0"], color="#0072B2", lw=2.0, label="MRG / EIG (n=0)")

    # Rossby and inertia-gravity modes
    colors = ["#3B5778", "#64809D", "#9AAFC1"]

    for i, n in enumerate(meridional_modes):
        color = colors[i % len(colors)]

        ax.plot(wn, curves[f"ER_{n}"], color=color, lw=1.5, label=f"ER (n={n})")
        ax.plot(wn, curves[f"EIG_{n}"], color=color, lw=1.5, ls="--", label=f"EIG (n={n})")
        ax.plot(wn, curves[f"WIG_{n}"], color=color, lw=1.5, ls=":", label=f"WIG (n={n})")

    # Reference axes
    ax.axvline(0, color="0.5", lw=0.8)
    ax.axhline(0, color="0.5", lw=0.8)

    # Axes limits and labels
    ax.set_xlim(-max_wn, max_wn)
    ax.set_ylim(0, max_freq)
    ax.set_xlabel("Zonal wavenumber", fontsize=text_size)
    ax.set_ylabel("Frequency (cycles day$^{-1}$)", fontsize=text_size)

    # Title
    ax.set_title(f"Matsuno dispersion curves  |  Equivalent depth = {he:g} m", fontsize=text_size+2, pad=15)

    # Figure styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", labelsize=text_size-1)

    # Legend in lower-right corner
    ax.legend(loc="lower right", fontsize=text_size-3, ncol=2, frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)

    # Save figure
    fig.tight_layout(rect=(0, 0.035, 1, 1))

    # Outside the axes to avoid covering the lower-right legend
    fig.text(0.975, 0.012, "PyWKSpectra", ha="right", va="bottom", fontsize=max(8, text_size-3), color="#243449", fontweight="bold")

    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()

    return fig, ax
