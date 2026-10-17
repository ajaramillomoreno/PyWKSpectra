# NOAA OLR Example Dataset

This directory contains a small, real-world dataset for demonstrating and testing **PyWKSpectra**, a Python package for Wheeler–Kiladis space–time spectral analysis.

## Included file

**`noaa_olr_daily_2000_2004.nc`** is a spatial and temporal subset of NOAA's daily interpolated outgoing longwave radiation (OLR) dataset.

| Property | Description |
| --- | --- |
| Variable | `olr` (outgoing longwave radiation, W m⁻²) |
| Period | 2000-01-01 to 2004-12-31 (1,827 daily time steps) |
| Latitude | 15°S to 15°N (13 points at 2.5° resolution) |
| Longitude | 0° to 357.5°E (144 points at 2.5° resolution; global zonal coverage) |
| Dimensions | `time: 1827`, `lat: 13`, `lon: 144` |
| Format | Compressed NetCDF4, approximately 5.30 MiB |
| Original NOAA file | `olr.day.mean.nc` |

The subset was created by selecting the specified dates and latitudes, retaining all longitudes and the original OLR values. **No anomalies, detrending, spectral filtering, regridding, or interpolation were applied during subset creation.** NetCDF compression is lossless.

The saved OLR values, coordinates, and dimensions were checked against the selected source data and matched exactly. This is a **demonstration and software-testing dataset**, not a claim of climatological representativeness or a complete validation of the Wheeler–Kiladis method.

## Data source and citation

Data provided by the **NOAA Physical Sciences Laboratory (PSL)**, Boulder, Colorado, USA.

- Dataset: [NOAA Interpolated Outgoing Longwave Radiation](https://psl.noaa.gov/data/gridded/data.interp_OLR.html)
- Provider: [NOAA PSL](https://psl.noaa.gov/)
- Original source filename: `olr.day.mean.nc`

**Recommended dataset reference:**

Liebmann, B., and C. A. Smith (1996): Description of a Complete (Interpolated) Outgoing Longwave Radiation Dataset. *Bulletin of the American Meteorological Society*, **77**(6), 1275–1277.

Please acknowledge NOAA PSL and cite the dataset reference when using these data in scientific publications. If your work also uses **PyWKSpectra**, cite the software separately according to the instructions in the project's main `README.md` and (when available) `CITATION.cff`.

## Reuse and disclaimer

NOAA PSL states that information on its government servers is generally in the **public domain unless otherwise noted**. This data subset is provided with attribution to NOAA; it is **not covered by the PyWKSpectra software license**. Neither this subset nor PyWKSpectra is an official NOAA product, and no NOAA endorsement or affiliation is implied.

See [NOAA PSL's usage and citation guidance](https://psl.noaa.gov/data/obs/download/Citation.html).
