"""Run parameters shared by the three notebooks.

Edit the values here, or override any of them with the environment variables
named below. The notebooks read and write each other's CSVs, so they must agree
on country, period and paths.
"""

from __future__ import annotations

import os
from pathlib import Path

import pycountry

# ── Study definition ──────────────────────────────────────────────────────────
ISO3 = os.getenv("TS_ISO3", "LCA")                       # any ISO-3166 alpha-3
BASIN = os.getenv("TS_BASIN", "NA")                      # IBTrACS basin code
START_YEAR = int(os.getenv("TS_START_YEAR", 1940))
END_YEAR = int(os.getenv("TS_END_YEAR", 2025))           # last complete season
IMPACT_RADIUS_KM = float(os.getenv("TS_RADIUS_KM", 150))

_country = pycountry.countries.get(alpha_3=ISO3)
if _country is None:
    raise ValueError(f"TS_ISO3={ISO3!r} is not a valid ISO-3166 alpha-3 code")
COUNTRY_NAME = _country.name
RECORD_YEARS = END_YEAR - START_YEAR + 1

# ── Impact criteria (knots — IBTrACS native unit) ─────────────────────────────
HURRICANE_MIN_KT = 64.0        # Saffir-Simpson Cat 1
TS_MIN_KT = 50.0               # severe TS floor; use 34.0 for all TS-force winds
TS_RADIUS_REDUCTION_KM = 50.0  # tighter buffer for sub-hurricane systems

IMPACT_CRITERIA = dict(
    hurricane_min_kt=HURRICANE_MIN_KT,
    ts_min_kt=TS_MIN_KT,
    radius_reduction_km=TS_RADIUS_REDUCTION_KM,
)

# ── Rainfall extraction (ERA5-Land via Earth Engine) ──────────────────────────
GEE_PROJECT = os.getenv("EE_PROJECT", "")    # required by notebook 02
ERA5_COLLECTION = "ECMWF/ERA5_LAND/HOURLY"
ERA5_BAND = "total_precipitation_hourly"     # hourly rate in m; ×1000 → mm/h
ERA5_SCALE_M = 11132                         # ~0.1° at the equator
GEOMETRY_BUFFER_DEG = 0.05                   # ~5 km, so small islands catch a cell

# Window: total hours centred on closest approach, and whether to also compute
# grid-max total and peak hourly intensity.
WINDOW_CONFIGS = {
    "3d": (72, False),
    "5d": (120, False),
    "7d": (168, True),
}

# ── Candidate AA thresholds (indicative — set per country) ────────────────────
THRESH_3D_MM = 100
THRESH_5D_MM = 150
THRESH_PEAK_MMH = 10

# ── EM-DAT ────────────────────────────────────────────────────────────────────
EMDAT_EXCLUDE_TYPES = ["Epidemic", "Mass movement (wet)"]   # not weather-driven
EMDAT_MATCH_TOLERANCE_DAYS = 5                              # EM-DAT ↔ IBTrACS

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / os.getenv("TS_INPUT_DIR", "inputs")
OUTPUT_DIR = ROOT / os.getenv("TS_OUTPUT_DIR", "outputs")

IBTRACS_CSV = INPUT_DIR / f"ibtracs.{BASIN}.list.v04r01.csv"
EMDAT_XLSX = INPUT_DIR / f"emdat_{ISO3.lower()}.xlsx"
MET_EVENTS_CSV = INPUT_DIR / f"{ISO3.lower()}_met_rainfall_events.csv"   # optional

STORMS_CSV = OUTPUT_DIR / f"ibtracs_impact_storms_{ISO3}_{START_YEAR}_{END_YEAR}.csv"
RAINFALL_CSV = OUTPUT_DIR / f"{ISO3}_rainfall_era5_{START_YEAR}_{END_YEAR}.csv"


def ensure_dirs():
    """Create the input/output directories if missing."""
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def summary():
    """Parameter block to print at the top of a notebook."""
    return "\n".join([
        f"Country  : {COUNTRY_NAME} ({ISO3})",
        f"Basin    : {BASIN}",
        f"Period   : {START_YEAR}–{END_YEAR}  ({RECORD_YEARS} years)",
        f"Radius   : {IMPACT_RADIUS_KM:.0f} km",
        f"Inputs   : {INPUT_DIR}",
        f"Outputs  : {OUTPUT_DIR}",
    ])
