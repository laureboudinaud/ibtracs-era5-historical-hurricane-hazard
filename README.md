# Historical Tropical Cyclone Hazard Analysis — Small Island States

Reproducible analysis of historical tropical-cyclone hazard and impact for a
single country, built from three open datasets:

| Pipeline | Source | Question it answers |
|---|---|---|
| **Wind & tracks** | IBTrACS v04r01 | Which storms passed close enough, and hard enough, to matter? How often? |
| **Rainfall** | ERA5-Land hourly (via Google Earth Engine) | How much rain did those storms deliver, over 3, 5 and 7 days? |
| **Impact** | EM-DAT | What did those events cost in people affected, deaths and damage? |

The outputs are intended as evidence for **anticipatory action (AA) trigger
design**: empirical exceedance curves and return periods a threshold can be
anchored to, rather than a threshold picked by convention.

Written for Eastern Caribbean small island developing states, but the code is
country-agnostic and driven by an ISO3 code.

---

## Repository layout

```
.
├── config.py                        # country, period, thresholds, paths
├── storm_utils.py                   # IBTrACS loading, geometry, distance, per-storm metrics
├── plot_utils.py                    # shared style and reusable figures
├── conftest.py                      # puts the repo root on sys.path for pytest
├── notebooks/
│   ├── 01_windspeed_ibtracs.ipynb   # run first  → outputs/ibtracs_impact_storms_*.csv
│   ├── 02_rainfall_era5.ipynb       # run second → outputs/*_rainfall_era5_*.csv
│   └── 03_impact_emdat.ipynb        # run third  → cross-references the two above
├── tests/test_utils.py
├── inputs/                          # data you download yourself (git-ignored)
├── outputs/                         # generated CSVs and figures (git-ignored)
└── environment.yml
```

The notebooks are **sequential**: `02` reads the storm list written by `01`, and
`03` joins EM-DAT against it. Out of order they fail on a missing CSV.

---

## Installation

```bash
conda env create -f environment.yml
conda activate ts-hazard
python -m pytest tests/ -q          # optional sanity check
```

`pip install -r requirements.txt` also works, but `cartopy` and `geopandas`
install more reliably from conda-forge.

The notebooks put the repo root on `sys.path` themselves, so they run from
either `notebooks/` or the repo root without installing anything.

Earth Engine (notebook `02` only):

```bash
earthengine authenticate
export EE_PROJECT=your-gee-project-id     # see .env.example
```

Notebook outputs are stripped before committing:

```bash
nbstripout --install                # honours the filter in .gitattributes
```

---

## Data

No input data is committed — download it into `inputs/` yourself.

**IBTrACS v04r01** — NOAA NCEI. Download the full CSV for your basin:
<https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/>

Use the NCEI file rather than the HDX mirror: the HDX extract carries only
`WMO_WIND`, which is missing for many storm-hours, and the fallback chain across
agency wind columns is what makes the record usable back to the mid-20th century.

**ERA5-Land hourly** — via Earth Engine (`ECMWF/ERA5_LAND/HOURLY`), no manual
download. Band `total_precipitation_hourly` is an hourly *rate* in metres, so it
needs no de-accumulation; ×1000 for mm h⁻¹.

**EM-DAT** — CRED / UCLouvain. Free account and a custom request at
<https://public.emdat.be/data>. **EM-DAT's terms do not permit redistribution**,
so the `.xlsx` must stay out of the repository; only derived aggregates are
shareable, with citation.

**Administrative boundaries** — geoBoundaries ADM0, fetched at runtime and cached
in `inputs/boundaries/` (CC BY 4.0).

**Met-service rainfall events (optional)** — if a national met service has shared
a list of non-cyclone rainfall events, save it as
`inputs/{iso3}_met_rainfall_events.csv` with at least `event_id`,
`event_datetime` and `name`, plus `ss_category` (use `Non-TC`) and `obs_rain_mm`
if available. Notebook `02` picks it up automatically and validates ERA5 against
`obs_rain_mm`; without the file it runs on IBTrACS alone.

---

## Configuration

Everything runtime-variable lives in `config.py` and can be overridden by
environment variable (`TS_ISO3`, `TS_BASIN`, `TS_START_YEAR`, `TS_END_YEAR`,
`TS_RADIUS_KM`, `TS_INPUT_DIR`, `TS_OUTPUT_DIR`, `EE_PROJECT`):

```python
ISO3              = "LCA"    # country
BASIN             = "NA"     # IBTrACS basin
START_YEAR        = 1940
END_YEAR          = 2025
IMPACT_RADIUS_KM  = 150      # proximity buffer around the coastline
HURRICANE_MIN_KT  = 64       # Saffir-Simpson Cat 1
TS_MIN_KT         = 50       # severe tropical storm floor
```

A storm enters the impact set when **either**:

* it is at hurricane force (≥ 64 kt) within `IMPACT_RADIUS_KM` of the coast, **or**
* it is at 50–63 kt within `IMPACT_RADIUS_KM − 50` km — a tighter buffer, on the
  basis that weaker systems have a smaller damaging-wind field.

Both are configurable; the 50 kt floor is a deliberate choice to exclude marginal
tropical storms, not the 34 kt TS definition.

---

## Outputs

| File | Contents |
|---|---|
| `ibtracs_impact_storms_{ISO3}_{start}_{end}.csv` | One row per qualifying storm: genesis position, closest approach, max wind, category, hours at hurricane/TS force, travel days |
| `{ISO3}_rainfall_era5_{start}_{end}.csv` | Per-event 3/5/7-day rainfall, 7-day grid maximum, peak hourly intensity |
| `{ISO3}_rainfall_era5_summary_statistics.csv` | Distribution summary across all events |
| `summary_by_hazard_{ISO3}.csv` | EM-DAT aggregates by hazard type |
| `outputs/*.png` | Exceedance curves, distributions, rainfall-vs-category and rainfall-vs-distance scatters, formation-zone maps, decade aggregations |

---

## Method notes and known limitations

These constrain how far the results can be pushed.

**Wind speeds mix reporting agencies.** `WIND` takes the first non-null value
across the agency columns (WMO first, then USA, …). Agencies do not share a
wind-averaging period — WMO reports 10-minute sustained winds, the US agencies
1-minute — so the series is not homogeneous, and Saffir-Simpson is defined on
1-minute winds. `WIND_SOURCE` records which agency supplied each value, so this
can be audited or filtered.

**Categories are assigned in knots.** IBTrACS winds come in 5-knot bins.
Classifying on a converted km/h value pushes a 64 kt storm to 118.5 km/h, below
the 119 km/h Cat 1 line, and silently demotes it to TS. Thresholds are applied in
knots and converted only for display; `saffir_simpson_category()` requires an
explicit `units` argument for the same reason.

**Track spacing is not always 6-hourly.** IBTrACS v04 interleaves 3-hourly
records, so duration metrics (`hours_as_hurricane`, `hours_as_ts`) use the median
observed timestep of each track rather than an assumed 6 hours.

**Distances are planar.** `add_distance` measures in a projected CRS. Use
`storm_utils.local_aeqd_crs(gdf)` — azimuthal-equidistant, centred on the country
— rather than a fixed UTM zone, which distorts badly once storms are several
hundred kilometres from the central meridian.

**ERA5-Land is coarse for small islands.** At 0.1° (~9 km) and masked to land, a
small island is covered by a handful of pixels. The area mean smooths away the
orographic maxima that drive flash flooding and will systematically understate
the local extreme. Treat these values as an island-scale index, not a point
rainfall estimate, and cross-check against station data where any exists.
Notebook `02` prints the cell count over the geometry; at one or two cells the
grid-maximum column is not a second metric.

**Rainfall windows are centred on closest approach.** A "3-day total" spans 36
hours either side. That characterises the event's total water; it is *not* a
forecast-window accumulation, so an AA trigger calibrated on these numbers must
be re-expressed in forecast terms before operational use.

**Return periods are empirical.** The observed rate over the record length, with
no distribution fitted and no confidence interval. With ~20 events the upper
percentiles rest on two or three storms. Exceedance curves use the Weibull
plotting position, so the largest observed event is not assigned a zero
exceedance probability, but nothing here extrapolates beyond the record.

**EM-DAT is inconsistent through time.** Reporting completeness improves sharply
after ~1980, so decade comparisons partly measure reporting effort rather than
hazard. The cross-reference matches on start-date proximity within
`EMDAT_MATCH_TOLERANCE_DAYS` and falls back to a year-only match, flagged as
ambiguous, where EM-DAT gives no month or day — check those by hand.

**Trend lines are descriptive.** The sample is conditioned on storms that passed
the proximity and intensity filters, and the reanalysis is not homogeneous across
the pre-satellite era. A fitted slope here is not evidence of a climate trend.

---

## Citation

If you use this code or its outputs, please cite the underlying data:

* Gahtan, J., K. R. Knapp, C. J. Schreck, H. J. Diamond, J. P. Kossin, M. C. Kruk (2024).
  *International Best Track Archive for Climate Stewardship (IBTrACS) Project, Version 4r01.*
  NOAA National Centers for Environmental Information.
* Muñoz Sabater, J. (2019). *ERA5-Land hourly data.* Copernicus Climate Change
  Service (C3S) Climate Data Store. Contains modified Copernicus Climate Change
  Service information.
* EM-DAT, CRED / UCLouvain, Brussels, Belgium — <https://www.emdat.be>
* Runfola, D. et al. (2020). *geoBoundaries: A global database of political
  administrative boundaries.* PLoS ONE 15(4).

---

## Licence

Code released under the MIT Licence (see `LICENSE`). Input datasets remain under
their own terms — in particular EM-DAT data may not be redistributed.

## Disclaimer

This repository contains analysis code and derived statistics only. It is
provided as-is, carries no operational warranty, and does not represent the
position of any institution. Thresholds shown are illustrative and must be
validated locally before any operational use.
