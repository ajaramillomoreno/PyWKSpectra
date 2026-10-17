<p align="center">
  <img src="https://raw.githubusercontent.com/ajaramillomoreno/PyWKSpectra/main/assets/pywkspectra_logo.png"
       alt="PyWKSpectra logo"
       width="280">
</p>


# PyWKSpectra
**Wheeler–Kiladis Space–Time Spectral Analysis in Python**

*A focused Python package for investigating tropical atmospheric variability through wavenumber–frequency spectra.*

**Current version:** 1.0

## About PyWKSpectra

**PyWKSpectra** provides a practical, reproducible way to compute and visualize Wheeler–Kiladis space–time spectra, a widely used diagnostic for investigating convectively coupled equatorial waves and tropical atmospheric variability.

Its philosophy is simple: **make the scientific analysis accessible without hiding the results**. The package combines the main processing steps in an easy-to-use interface while returning a standard **xarray Dataset**, so researchers can inspect the spectral estimates, save them, and create their own analyses and figures.

PyWKSpectra 1.0 includes:

- Input handling for NetCDF, xarray objects, and GRIB data (with optional support), including coordinate and sampling checks.
- Removal of dominant low-frequency variability, equatorial symmetric/antisymmetric decomposition, and windowed space–time Fourier analysis.
- Symmetric and antisymmetric power spectra, smoothed background power, and background-normalized spectra.
- Ready-to-use plots of normalized and raw spectra, plus theoretical Matsuno dispersion diagrams.
- NetCDF-compatible results that can be reused without repeating the analysis.

**Current scope:** PyWKSpectra computes and visualizes spectral diagnostics. **Version 1.0 does not perform Wheeler–Kiladis wavenumber–frequency filtering or reconstruct filtered wave signals in physical space.** Those capabilities, along with other improvements, remain possible directions for future development.

## What is new compared with NCL and `wk_spectra`?

PyWKSpectra builds on the Wheeler–Kiladis methodology and on the experience gained from NCL's `wkSpaceTime` and my earlier Python implementation, `wk_spectra`. **It is not a line-by-line translation of NCL or simply a renamed version of the previous package.** Version 1.0 adopts a more Python-native approach, with an emphasis on clarity, reproducibility, and a manageable set of direct dependencies.

Some of the main changes and improvements are:

- **An xarray-centered workflow.** PyWKSpectra accepts NetCDF files and xarray `Dataset` or `DataArray` objects, with optional GRIB support. It returns a labeled `xarray.Dataset` containing spectral coordinates, eight output variables, and analysis metadata. Results can be explored directly or saved as NetCDF for later use.
- **Simpler, more transparent interfaces.** One public function, `compute_wk_spectra()`, coordinates the complete analysis. Data input, spectral computations, and plotting are organized into separate modules rather than requiring users to manage internal plotting state.
- **Explicit data checks and preprocessing.** The package checks regular temporal sampling, complete global longitude coverage, and symmetric latitude pairs. It documents its detrending, low-frequency removal, overlapping windows, and spectral-smoothing choices.
- **NumPy-based numerical calculations.** Space–time Fourier analysis uses NumPy FFT routines, and preprocessing can operate in latitude chunks to manage working memory. The implementation uses a Hann window for each spectral segment; this is a methodological choice and need not produce numerically identical results to NCL or older versions of `wk_spectra`.
- **Integrated scientific figures.** Standard functions generate raw spectra, a two-panel background-normalized spectrum, and Matsuno dispersion diagrams. The dispersion curves are calculated with NumPy, without requiring SciPy root-finding or pandas DataFrames in the plotting code.
- **Fewer direct dependencies.** The core package specifies NumPy, xarray, Matplotlib, and netCDF4. SciPy and pandas are not direct requirements of PyWKSpectra's own algorithms, although **xarray may install pandas as a dependency**. GRIB support remains optional through `cfgrib` and ecCodes.
- **A reproducible starting point.** A documented notebook and a small NOAA OLR dataset make it easier to learn the workflow, inspect outputs, and adapt the analysis to other datasets.

These improvements concern **software design, usability, and reproducibility**; they are not a claim that PyWKSpectra is universally faster, more accurate, or numerically equivalent to the NCL diagnostic or earlier `wk_spectra` releases. Scientific comparisons will continue as the package evolves. Likewise, PyWKSpectra 1.0 is focused on spectral diagnostics, not wavenumber–frequency filtering or reconstruction of individual waves.

## Installation and examples

### Install the package

PyWKSpectra requires **Python 3.10 or newer**. Install it with pip:

```bash
pip install PyWKSpectra
```

The core scientific dependencies are installed automatically:

| Package | Purpose |
|---|---|
| **NumPy** | Numerical operations and Fourier transforms |
| **xarray** | Multidimensional atmospheric data and spectral outputs |
| **Matplotlib** | Spectral and dispersion-diagram visualization |
| **netCDF4** | NetCDF input/output support |

**Optional GRIB support.**  
GRIB decoding also requires **ecCodes**. Depending on your operating system and environment, its native libraries may need to be installed separately; for example, with Conda:

```bash
conda install -c conda-forge eccodes
```

### Run the example notebook

The repository includes a step-by-step tutorial:

[`examples/example_wk_spectra.ipynb`](examples/example_wk_spectra.ipynb)

It uses the accompanying **NOAA daily outgoing longwave radiation (OLR)** dataset for **2000–2004**:

[`examples/data/noaa_olr_daily_2000_2004.nc`](examples/data/noaa_olr_daily_2000_2004.nc)

To run the notebook, you will also need a Jupyter-compatible environment and **ipykernel**:

```bash
python -m pip install ipykernel
```

You can open the notebook using JupyterLab or Visual Studio Code with its Python and Jupyter extensions. Select the Python kernel in which PyWKSpectra and `ipykernel` are installed. Jupyter tools are needed for the tutorial, **not** for running PyWKSpectra in regular Python scripts.

### Quick start

```python
import pywkspectra as wk

results = wk.compute_wk_spectra(
    "examples/data/noaa_olr_daily_2000_2004.nc",
    varname="olr",
    latBound=15,
)

# Display symmetric and antisymmetric background-normalized spectra
fig, axes = wk.plot_wk_spectra(results)

# Save spectral results for later use
results.to_netcdf("wk_spectra_results.nc")
```

Additional plotting functions include `wk.plot_raw_spectra(results)` and `wk.plot_matsuno(he=25)`. To save figures in PNG and PDF, pass a directory using `save_dir`, for example `wk.plot_wk_spectra(results, save_dir="figures", dpi=300)`.

## Documentation

PyWKSpectra provides scientific documentation and a reproducible example to help researchers understand and apply Wheeler–Kiladis spectral analysis.

- **[Scientific methodology](docs/methodology.md)** — Detailed description of the spectral analysis procedure, including preprocessing, equatorial symmetry decomposition, Fourier transforms, spectral smoothing, background estimation, normalization, and theoretical dispersion relationships.

- **[Example notebook](examples/example_wk_spectra.ipynb)** — A step-by-step demonstration using real NOAA outgoing longwave radiation (OLR) data, including spectral calculations, visualization, and export of results.

The methodology document describes the scientific and numerical foundations of the package, while the example notebook provides a practical introduction to its application.

## A little history (and a note on attribution)

This project began long before the name **PyWKSpectra** existed.

I started working on a Python implementation of Wheeler–Kiladis spectral analysis in **2014**, during my doctoral studies, while visiting **George Kiladis** at NOAA's Physical Sciences Division (PSD). My initial goal was to bring to Python a method I had been using through the **NCAR Command Language (NCL)**, especially its [`wkSpaceTime` diagnostic](https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml).

The NCL implementation was an invaluable scientific and technical reference. Over time, however, the Python code took its own path: it evolved toward NumPy-based calculations, a more Python-native workflow, and a package known as **`wk_spectra`**. An example of its subsequent development and circulation can be found in this [repository from the `wk_spectra` lineage](https://github.com/ajaramillomoreno/wk_spectra).

Like many research-software projects, `wk_spectra` was developed alongside demanding academic responsibilities and without dedicated time or funding for software maintenance. As interest in the code grew, I could not always provide the documentation, user support, and updates it deserved. I was also disappointed to see that some later uses and citations did not correctly identify me as the original developer, despite the substantial work that went into the implementation.

I mention this not to criticize individual researchers or overlook later contributions, but because **software provenance and proper attribution matter just as much as scientific reproducibility**. Code can be copied, modified, and improved by many people; its development history should remain clear.

**PyWKSpectra** is a fresh effort to modernize and clarify the software, provide better documentation, establish its authorship explicitly, and create a foundation for future scientific development. Contributions by others to earlier versions deserve appropriate acknowledgment, too.

The project continues to be developed **in my limited personal time, without dedicated funding or institutional support for its development**, alongside substantial research, teaching, academic, and administrative responsibilities. My institutional affiliation identifies where I work, not a funding source for this project. I welcome questions, suggestions, and reasonable requests, but please understand that replies may be delayed—or occasionally not possible. Thank you for your patience.

## How to cite PyWKSpectra

**If PyWKSpectra contributes to your figures, analyses, teaching materials, or scientific results, please cite the software and its original author.**

PyWKSpectra 1.0 is officially archived on **Zenodo**, providing a permanent DOI for scientific citation.

**Recommended citation:**

Jaramillo-Moreno, A. (2026). *PyWKSpectra: Wheeler–Kiladis Space–Time Spectral Analysis in Python* (Version 1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23438592

**Software DOI:** [10.5281/zenodo.23438592](https://doi.org/10.5281/zenodo.23438592)

The repository also includes a [CITATION.cff](CITATION.cff) file with machine-readable citation metadata.

The foundational paper by **Wheeler and Kiladis (1999)** describes the scientific methodology. A separate software citation recognizes the development, implementation, documentation, and maintenance that make this analysis available in Python. Ideally, please cite both when relevant.

## Publications using PyWKSpectra

**Have you used PyWKSpectra in a paper, preprint, thesis, or report? Please send me the reference or DOI!**

I would be delighted to learn about your work and maintain a list of publications that have benefited from the package. A publication list will be added as references become available.

## The best way to say thanks

I am not asking for donations or payment. If PyWKSpectra has helped your research, three things would mean a great deal to me:

1. **Cite the software and acknowledge its original development.**
2. **Share your publication** so I can learn from it and add it to the list of work using PyWKSpectra.
3. **Invite me to collaborate** when there is a meaningful opportunity to explore tropical atmospheric dynamics together.

Bug reports, constructive feedback, and properly credited contributions are welcome as well. And if you use the package in a classroom, for a student project, or simply for the joy of exploring atmospheric waves, I would be just as happy to know that these hours of work have been useful to someone.

*Coffee is appreciated, but citations and scientific collaborations are even better fuel for research software.*

## Scientific references

- **Matsuno, T. (1966).** Quasi-geostrophic motions in the equatorial area. *Journal of the Meteorological Society of Japan*, **44**, 25–43. [https://doi.org/10.2151/jmsj1965.44.1_25](https://doi.org/10.2151/jmsj1965.44.1_25)
- **Wheeler, M., & Kiladis, G. N. (1999).** Convectively coupled equatorial waves: Analysis of clouds and temperature in the wavenumber–frequency domain. *Journal of the Atmospheric Sciences*, **56**, 374–399. [https://doi.org/10.1175/1520-0469(1999)056%3C0374:CCEWAO%3E2.0.CO;2](https://doi.org/10.1175/1520-0469(1999)056%3C0374:CCEWAO%3E2.0.CO;2)
- **Kiladis, G. N., Wheeler, M. C., Haertel, P. T., Straub, K. H., & Roundy, P. E. (2009).** Convectively coupled equatorial waves. *Reviews of Geophysics*, **47**, RG2003. [https://doi.org/10.1029/2008RG000266](https://doi.org/10.1029/2008RG000266)
- **NCL reference:** [`wkSpaceTime`](https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml).

## Author and maintainer

**Alejandro Jaramillo-Moreno**  
Departamento de Cambio Climático  
Instituto de Ciencias de la Atmósfera y Cambio Climático (ICAyCC)  
Universidad Nacional Autónoma de México (UNAM)  
Mexico City, Mexico

**ORCID:** [0000-0002-4175-7726](https://orcid.org/0000-0002-4175-7726)  
**Email:** [ajaramillo@atmosfera.unam.mx](mailto:ajaramillo@atmosfera.unam.mx)

## License and acknowledgments

**PyWKSpectra 1.0 is distributed under the GNU General Public License version 3 (GPL-3.0-only).** See the [LICENSE](LICENSE) file for the complete terms.

Copyright © 2026 Alejandro Jaramillo-Moreno.

PyWKSpectra builds upon the scientific and software development work originally undertaken by Alejandro Jaramillo, including the earlier `wk_spectra` implementation.

The original `wk_spectra` code was distributed under the BSD 2-Clause License, with the copyright notice:

*Copyright (c) 2018, Alejandro Jaramillo.*

Any portions of that code incorporated into PyWKSpectra retain their applicable original license terms and copyright notices. A copy of the original license is provided in [LICENSES/BSD-2-Clause-wk_spectra.txt](LICENSES/BSD-2-Clause-wk_spectra.txt).

PyWKSpectra also acknowledges the foundational scientific contributions of Wheeler and Kiladis (1999), the NCL `wkSpaceTime` implementation that inspired the original Python work, and the researchers who have contributed ideas, feedback, and improvements over the years.

The project is committed to open scientific software, transparent development, reproducibility, and proper recognition of scientific contributions.

---

*May your spectral peaks be physically meaningful.*
