# PyWKSpectra 1.0 — Scientific Methodology

**Wheeler–Kiladis space–time spectral analysis in Python**  
**Author:** Alejandro Jaramillo-Moreno  
Departamento de Cambio Climático, Instituto de Ciencias de la Atmósfera y Cambio Climático (ICAyCC), Universidad Nacional Autónoma de México (UNAM)

## 1. Purpose and scope

PyWKSpectra computes wavenumber–frequency power spectra for tropical atmospheric fields, following the general diagnostic framework of Wheeler and Kiladis (1999). It separates variability by symmetry about the equator, estimates a smoothed spectral background, and provides both unnormalized and background-normalized spectra. It can also display the theoretical equatorial-wave dispersion relationships of Matsuno (1966).

The package is intended for analysis of variables such as outgoing longwave radiation (OLR), precipitation, or other regularly sampled tropical fields on a global longitude grid. **Version 1.0 does not perform Wheeler–Kiladis spectral filtering, inverse transforms, or reconstruction of individual waves in physical space.**

The description below documents the current numerical choices of PyWKSpectra 1.0. It does not imply exact numerical equivalence to NCL `wkSpaceTime`, the previous `wk_spectra` implementation, or published spectra produced with other parameters.

## 2. Input data and coordinate handling

The main public function is:

```python
import pywkspectra as wk

results = wk.compute_wk_spectra(source, varname="olr", latBound=15)
```

`source` may be a NetCDF or GRIB path, or an `xarray.Dataset`/`xarray.DataArray`. GRIB reading requires the optional `cfgrib`/ecCodes backend. The input variable must be a three-dimensional field ordered logically by **time, latitude, and longitude**. Common coordinate names are recognized and standardized by `load_wk_data()`.

The input preparation requires:

- Strictly increasing, equally spaced timestamps, with a sampling interval that divides one day exactly. The sampling frequency, `spd`, is expressed as samples per day.
- A regular, complete 360° longitude grid. Longitudes are standardized to `[-180°, 180°)` and sorted in increasing order.
- Latitude pairs symmetric about the equator. The requested domain is `-latBound` to `+latBound`; if the available coverage is asymmetric, the function reduces the bound to the largest available symmetric subset. Latitude is ordered from north to south.
- Real, finite values over the selected domain. Missing and non-finite values are not interpolated or silently filled.

The selected data are loaded into memory in version 1.0. Users should consider memory requirements for long records and high-resolution grids.

## 3. Time windows and sampling

Let `N_t` denote the number of input time samples, `s` the samples per day, `T` the window length (`nDayWin`, in days), and `D` the window skip (`nDaySkip`, in days). The window length and advance in samples are

$$
N = sT, \qquad A = s(T+D).
$$

A **negative** `nDaySkip` gives overlap; for example, `nDayWin=96` and `nDaySkip=-48` produce 96-day windows advancing every 48 days, with 48-day overlap. Provided there are sufficient samples, the number of complete windows is

$$
N_{\mathrm{win}} = 1 + \left\lfloor \frac{N_t-N}{A} \right\rfloor.
$$

Any incomplete remainder at the end of the record is not used. For regularly sampled data, the nominal frequency resolution within each window is

$$
\Delta f = \frac{s}{N} = \frac{1}{T}\quad \text{cycles day}^{-1}.
$$

For daily observations and 96-day windows, `Δf = 1/96 ≈ 0.01042` cycles per day. The Nyquist frequency is `s/2` cycles per day.

## 4. Removal of dominant temporal signals

Before splitting the record into windows, PyWKSpectra applies temporal preprocessing independently at each grid point:

1. Fit and remove a linear trend over the full input record, **preserving the grid-point temporal mean**.
2. For records of at least 365 days, apply a Fourier high-pass operation to the full time series: nonzero temporal-frequency components with periods **greater than or equal to the spectral window length** are set to zero. The zero-frequency mean is retained.

With a 96-day window, the nominal cutoff is

$$
f_c = \frac{1}{96}\quad \text{cycles day}^{-1}.
$$

This is a Fourier-domain cutoff applied to the complete record, not a separately estimated seasonal-cycle climatology. The sharp frequency cutoff and treatment of the record endpoints may influence behavior near the cutoff; this should be considered when interpreting the lowest retained frequencies. For shorter records, the full-record Fourier high-pass stage is not applied, although detrending is still performed.

## 5. Equatorial symmetry decomposition

At each time, longitude, and latitude magnitude $\phi>0$, let $X_N=X(+\phi)$ and $X_S=X(-\phi)$. The symmetric and antisymmetric components are

$$
X_{\mathrm{sym}}(\phi)=\frac{X_N+X_S}{2},
\qquad
X_{\mathrm{anti}}(\phi)=\frac{X_N-X_S}{2}.
$$

The implementation stores the **antisymmetric** component at the positive (northern) latitudes and the **symmetric** component at the negative (southern) latitudes within one array. This is an internal storage convention, not a geographic restriction on either type of wave.

The equator, when present, is retained in the intermediate data but **excluded from the latitude average** used to form the final spectra. The subsequent spectral averaging uses corresponding latitude pairs in the selected tropical domain.

The interpretation of symmetry depends on the analyzed physical variable; for example, conventional symmetry labels for scalar OLR need not be applied uncritically to every wind component.

## 6. Space–time Fourier analysis

Each overlapping window is processed separately. A linear temporal trend and temporal mean are removed within the window, and a **Hann taper** is applied along its time axis:

$$
w_j=\tfrac12\left[1-\cos\left(\frac{2\pi j}{N-1}\right)\right],
\qquad j=0,\ldots,N-1.
$$

A two-dimensional discrete Fourier transform is then calculated over **time and longitude**. In symbolic form, the windowed transform is

$$
\widehat X(f,k)=\frac{1}{N N_\lambda}
\sum_{j=0}^{N-1}\sum_{\ell=0}^{N_\lambda-1}
\left[X_{j,\ell}-\widehat{\text{linear trend}}_{j,\ell}\right]
\,w_j\,\exp\!\left[-2\pi i\left(\frac{mj}{N}+\frac{n\ell}{N_\lambda}\right)\right],
$$

where $N_\lambda$ is the number of longitudes and $m,n$ identify discrete Fourier bins. The code uses `numpy.fft`, rearranges the transform into its analysis wavenumber–frequency convention, and retains positive temporal frequencies in the exported result.

The **power** at each spectral bin is computed as

$$
P(f,k)=\left|\widehat X(f,k)\right|^2.
$$

The normalization above describes the transform scaling used in the implementation. The output is **squared Fourier amplitude**, not a separately calibrated power spectral density per unit frequency. In particular, the Hann taper is **not** subsequently corrected by an explicit taper-energy factor. Consequently, absolute power magnitudes should only be compared across analyses with attention to their sampling, taper, and normalization settings.

The plotting convention represents eastward-propagating features at **positive** zonal wavenumbers and westward-propagating features at **negative** wavenumbers. Because spectral direction depends on FFT sign and axis-reordering conventions, this convention should be verified with known propagating signals before relying on it for new data configurations.

## 7. Symmetric and antisymmetric power spectra

For each symmetry class, the code first computes the squared modulus of the Fourier coefficients, then averages spectral power across all complete time windows and the latitude points assigned to that symmetry class:

$$
\overline P_{\mathrm{sym}}(f,k)
=\left\langle P_{\mathrm{sym}}^{(r,\phi)}(f,k)\right\rangle_{r,\phi},
\qquad
\overline P_{\mathrm{anti}}(f,k)
=\left\langle P_{\mathrm{anti}}^{(r,\phi)}(f,k)\right\rangle_{r,\phi}.
$$

Here $r$ indexes time windows. The present implementation takes an **unweighted arithmetic mean** over the selected latitudes and windows, rather than applying area weights.

Both spectra receive **one pass** of a 1–2–1 smoother along the **positive-frequency axis**, for zonal wavenumbers $|k|\leq 27$. For interior points, the operator is

$$
S_i=\tfrac14 P_{i-1}+\tfrac12 P_i+\tfrac14 P_{i+1},
$$

with adapted endpoint weights. The zero-frequency bin is excluded from exported spectra.

## 8. Background spectrum

A background spectrum is estimated from power across both hemispheres (excluding the equator), averaged over latitude and spectral windows. It is smoothed more strongly than the symmetric and antisymmetric spectra to represent the broad background against which localized spectral peaks are evaluated.

The implemented smoothing sequence is:

1. **Along zonal wavenumber**, apply the 1–2–1 smoother within $|k|\leq 27$, with a frequency-dependent number of passes:

   | Frequency (cycles/day) | Passes |
   |---|---:|
   | $0<f<0.1$ | 5 |
   | $0.1\leq f<0.2$ | 10 |
   | $0.2\leq f<0.3$ | 20 |
   | $f\geq 0.3$ | 40 |

2. **Along positive temporal frequency**, apply 10 passes of the 1–2–1 smoother for each wavenumber within $|k|\leq 27$, using the available positive-frequency range up to 0.8 cycles/day.

This is a **smoothed empirical background**, not a theoretical red-noise model or an independent estimate of statistical significance.

## 9. Background normalization and exported products

The symmetric and antisymmetric spectra are normalized by the same background:

$$
P_{\mathrm{sym,norm}}(f,k)=\frac{\overline P_{\mathrm{sym}}(f,k)}{P_{\mathrm{bg}}(f,k)},
\qquad
P_{\mathrm{anti,norm}}(f,k)=\frac{\overline P_{\mathrm{anti}}(f,k)}{P_{\mathrm{bg}}(f,k)}.
$$

Values greater than one indicate enhancement relative to this background estimate; **they are not, by themselves, formal significance levels**. PyWKSpectra also calculates $\log_{10}$ of the unnormalized symmetric, antisymmetric, and background power (with nonpositive values represented as missing where applicable).

The main function returns an `xarray.Dataset` with dimensions `frequency` (cycles/day) and `wavenumber`, containing eight variables:

| Variable | Meaning |
|---|---|
| `symmetric_power` | Symmetric spectral power |
| `antisymmetric_power` | Antisymmetric spectral power |
| `background_power` | Smoothed background power |
| `symmetric_log_power` | $\log_{10}$ of symmetric power |
| `antisymmetric_log_power` | $\log_{10}$ of antisymmetric power |
| `background_log_power` | $\log_{10}$ of background power |
| `symmetric_normalized` | Symmetric power / background |
| `antisymmetric_normalized` | Antisymmetric power / background |

Exported data are restricted to **positive frequencies**, $f\leq$ `max_freq`, and $|k|\leq$ `max_wn`. The default limits are `max_freq=0.5` cycles/day and `max_wn=20`. Dataset attributes record information about the input variable, period, latitude bound, sampling and window parameters, taper, and processing choices.

The resulting dataset can be saved for later analyses:

```python
results.to_netcdf("wk_spectra_results.nc")
```

## 10. Theoretical Matsuno dispersion curves

PyWKSpectra can compute and plot the theoretical dispersion relationships of equatorial waves using an equatorial beta-plane approximation. With equivalent depth $h_e$, gravity $g$, Earth's radius $a$, and rotation rate $\Omega$,

$$
c=\sqrt{gh_e}, \qquad \beta=\frac{2\Omega}{a}, \qquad k=\frac{n_\lambda}{a},
$$

where $n_\lambda$ is the signed global zonal wavenumber used in the plots. Angular frequency is converted to cycles per day by multiplication by $86400/(2\pi)$.

The eastward Kelvin branch uses $\omega=ck$ for $k\geq 0$. The continuous mixed Rossby–gravity / $n=0$ eastward inertia–gravity branch is

$$
\omega_{\mathrm{MRG/EIG0}}=
\frac{ck+\sqrt{c^2k^2+4\beta c}}{2}.
$$

For meridional mode numbers $n\geq 1$, the equatorial shallow-water dispersion relation is

$$
\omega^3-\left[c^2k^2+(2n+1)\beta c\right]\omega-\beta c^2k=0.
$$

The implementation evaluates the appropriate real roots using NumPy expressions, without SciPy root-finding or pandas tables, to display equatorial Rossby and eastward/westward inertia–gravity branches. The dispersion curves are interpretive references; their presence on a power-spectrum plot does not constitute mode identification or a significance test.

## 11. Reproducibility, differences, and limitations

- **Compared with NCL:** PyWKSpectra uses a NumPy/xarray workflow and a **Hann** temporal window, rather than assuming numerical equivalence with every tapering, preprocessing, and smoothing choice of NCL `wkSpaceTime`.
- **Compared with earlier `wk_spectra`:** it provides an independent data-loading interface, named coordinates and metadata, an `xarray.Dataset` output, a simplified public API, and consolidated plotting tools. The implementation also uses fewer *direct* numerical dependencies, although packages such as xarray may have their own dependencies.
- **No missing-data infilling or regridding:** users must prepare suitable regularly sampled and globally periodic input fields before analysis.
- **No filtering or reconstruction:** version 1.0 does not isolate selected spectral regions and transform them back into physical space.
- **No built-in significance test:** normalized peaks indicate relative power enhancement, not a reported statistical p-value.
- **Validation scope:** qualitative similarity to spectra in the scientific literature is encouraging, but it is not a substitute for independent numerical benchmarks, directional tests with synthetic waves, or comprehensive automated regression tests. The initial reproducible example uses a daily 144-longitude grid; alternative sampling and longitude-grid configurations deserve additional checks.

A documented example with NOAA daily interpolated OLR for 2000–2004 is available at [`examples/example_wk_spectra.ipynb`](../examples/example_wk_spectra.ipynb). Its input dataset and provenance notes are in [`examples/data/`](../examples/data/).

## 12. References

- **Matsuno, T. (1966).** Quasi-geostrophic motions in the equatorial area. *Journal of the Meteorological Society of Japan*, **44**, 25–43. https://doi.org/10.2151/jmsj1965.44.1_25
- **Wheeler, M., and G. N. Kiladis (1999).** Convectively coupled equatorial waves: Analysis of clouds and temperature in the wavenumber–frequency domain. *Journal of the Atmospheric Sciences*, **56**, 374–399. https://doi.org/10.1175/1520-0469(1999)056%3C0374:CCEWAO%3E2.0.CO;2
- **Kiladis, G. N., M. C. Wheeler, P. T. Haertel, K. H. Straub, and P. E. Roundy (2009).** Convectively coupled equatorial waves. *Reviews of Geophysics*, **47**, RG2003. https://doi.org/10.1029/2008RG000266
- **NCL `wkSpaceTime` diagnostic:** https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml

---

**Software:** PyWKSpectra 1.0  
**Contact:** [ajaramillo@atmosfera.unam.mx](mailto:ajaramillo@atmosfera.unam.mx)  
**Citation:** See [`CITATION.cff`](../CITATION.cff) for the software's recommended citation.
