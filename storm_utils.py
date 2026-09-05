"""IBTrACS loading, country geometry, distance to shore, per-storm metrics.

Conventions: winds in knots (IBTrACS native) and converted to km/h only for
display; distances in km from a track point to the nearest point on the country
boundary (0 over land); timestamps are IBTrACS ISO_TIME (UTC).
"""

from __future__ import annotations

import io
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

__version__ = "0.3.0"

KT_TO_KMH = 1.852

#: Saffir-Simpson lower bounds in knots. The familiar km/h figures
#: (119/154/178/209/252) are rounded conversions of these.
SS_THRESHOLDS_KT = {5: 137.0, 4: 113.0, 3: 96.0, 2: 83.0, 1: 64.0}

HURRICANE_MIN_KT = 64.0     # Cat 1
DEFAULT_TS_MIN_KT = 50.0    # severe TS floor used by ts_mask


# ── UNITS AND CATEGORY ────────────────────────────────────────────────────────
def knots_to_kmh(knots):
    """Knots → km/h (scalars, Series or arrays)."""
    return knots * KT_TO_KMH


def kmh_to_knots(kmh):
    """km/h → knots."""
    return kmh / KT_TO_KMH


def saffir_simpson_category(intensity, units):
    """Saffir-Simpson category: 1–5, 'TS', or 'Unknown' if missing.

    ``units`` ('kt' or 'kmh') is required on purpose: classifying a converted
    km/h value silently demotes a 64 kt storm (118.5 km/h < 119) to TS.
    """
    if units not in ("kt", "kmh"):
        raise ValueError("units must be 'kt' or 'kmh'")
    if intensity is None or not np.isfinite(np.asarray(intensity, dtype=float)):
        return "Unknown"

    kt = float(intensity) if units == "kt" else kmh_to_knots(float(intensity))
    for cat in sorted(SS_THRESHOLDS_KT, reverse=True):
        if kt >= SS_THRESHOLDS_KT[cat]:
            return cat
    return "TS"


def saffir_simpson_label(intensity, units):
    """As above, but always a string: 'TS', 'Cat 1'…'Cat 5', 'Unknown'."""
    cat = saffir_simpson_category(intensity, units)
    return f"Cat {cat}" if isinstance(cat, int) else cat


# ── COUNTRY BOUNDARY ──────────────────────────────────────────────────────────
GEOBOUNDARIES_API = "https://www.geoboundaries.org/api/current/gbOpen/{iso3}/ADM0/"


def load_boundary(iso3, crs_meters=None, cache_dir="inputs/boundaries", timeout=60):
    """Download (and cache) an ADM0 boundary from geoBoundaries.

    Returns ``(GeoDataFrame in EPSG:4326, dissolved boundary in crs_meters)``;
    the second element is None unless ``crs_meters`` is given. Pass
    ``cache_dir=None`` to disable caching.
    """
    cache_path = Path(cache_dir) / f"{iso3}_ADM0.geojson" if cache_dir else None

    if cache_path and cache_path.exists():
        gdf = gpd.read_file(cache_path)
    else:
        meta = requests.get(GEOBOUNDARIES_API.format(iso3=iso3), timeout=timeout)
        meta.raise_for_status()
        payload = requests.get(meta.json()["gjDownloadURL"], timeout=timeout)
        payload.raise_for_status()
        gdf = gpd.read_file(io.BytesIO(payload.content))
        if cache_path:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            gdf.to_file(cache_path, driver="GeoJSON")

    gdf = gdf.to_crs("EPSG:4326")

    boundary = None
    if crs_meters is not None:
        # union_all(), not .iloc[0]: archipelagic states span several features.
        boundary = gdf.to_crs(crs_meters).union_all()

    return gdf, boundary


def local_aeqd_crs(gdf):
    """Azimuthal-equidistant CRS centred on ``gdf``.

    Distances are then accurate in every direction from the country centre,
    which a fixed UTM zone is not several hundred km off its meridian.
    """
    centroid = gdf.to_crs("EPSG:4326").union_all().centroid
    return f"+proj=aeqd +lat_0={centroid.y} +lon_0={centroid.x} +units=m +datum=WGS84 +no_defs"


# ── IBTRACS LOADER ────────────────────────────────────────────────────────────
#: Agency wind columns, in fallback priority order (first non-null wins).
WIND_COLS = [
    "WMO_WIND", "USA_WIND", "TOK_WIND", "CMA_WIND", "HKO_WIND",
    "KMA_WIND", "NEW_WIND", "REU_WIND", "BOM_WIND", "NAD_WIND",
    "WEL_WIND", "DS8_WIND", "TD6_WIND", "TD5_WIND", "NEU_WIND",
]


def load_ibtracs(file_path, basin="NA", start_year=1940, end_year=2025,
                 whole_track=True, wind_cols=WIND_COLS):
    """Load and clean an IBTrACS v04 CSV.

    Adds ``LAT``/``LON`` (numeric), ``DATE`` (UTC), ``YEAR``, ``WIND`` (knots)
    and ``WIND_SOURCE`` (the agency column that supplied ``WIND``).

    ``whole_track=True`` keeps every point of any storm with at least one point
    in ``basin``; False truncates basin-crossing storms and so corrupts
    genesis-derived metrics. ``WIND`` mixes agencies, which do not share a
    wind-averaging period (WMO 10-minute, US 1-minute) — hence ``WIND_SOURCE``.
    """
    df = pd.read_csv(file_path, skiprows=[1], low_memory=False,
                     na_values=[" "], keep_default_na=False)

    df["LAT"] = pd.to_numeric(df["LAT"], errors="coerce")
    df["LON"] = pd.to_numeric(df["LON"], errors="coerce")
    df["DATE"] = pd.to_datetime(df["ISO_TIME"], errors="coerce")
    df["YEAR"] = df["DATE"].dt.year
    df["SEASON"] = (pd.to_numeric(df["SEASON"], errors="coerce")
                    if "SEASON" in df.columns else np.nan)

    if whole_track:
        basin_sids = df.loc[df["BASIN"] == basin, "SID"].unique()
        df = df[df["SID"].isin(basin_sids)]
    else:
        df = df[df["BASIN"] == basin]

    # Filter on season, not point-year: storms can straddle 31 December.
    season = df["SEASON"].fillna(df["YEAR"])
    df = df[(season >= start_year) & (season <= end_year)]
    df = df.dropna(subset=["LAT", "LON", "DATE"]).copy()

    available = [c for c in wind_cols if c in df.columns]
    if not available:
        raise ValueError(f"None of the expected wind columns found: {wind_cols}")

    df[available] = df[available].apply(pd.to_numeric, errors="coerce")
    filled = df[available].notna()
    df["WIND"] = df[available].bfill(axis=1).iloc[:, 0]
    df["WIND_SOURCE"] = filled.idxmax(axis=1).where(filled.any(axis=1))
    df = df.dropna(subset=["WIND"])

    return df.sort_values(["SID", "DATE"]).reset_index(drop=True)


def add_distance(df, boundary, crs_meters):
    """Add ``distance_km``: track point → country boundary (0 over land)."""
    df = df.copy()
    df["distance_km"] = (
        gpd.GeoSeries.from_xy(df["LON"], df["LAT"], crs=4326)
        .to_crs(crs_meters)
        .distance(boundary)
        / 1000
    )
    return df


# ── IMPACT CRITERIA (single source of truth) ──────────────────────────────────
def hurricane_mask(group, impact_radius_km, hurricane_min_kt=HURRICANE_MIN_KT):
    """Timesteps at hurricane force within ``impact_radius_km``."""
    return (group["WIND"] >= hurricane_min_kt) & (group["distance_km"] < impact_radius_km)


def ts_mask(group, impact_radius_km, ts_min_kt=DEFAULT_TS_MIN_KT,
            hurricane_min_kt=HURRICANE_MIN_KT, radius_reduction_km=50.0):
    """Timesteps at sub-hurricane TS force within a tightened radius.

    Weaker systems must pass closer, their damaging-wind field being smaller.
    """
    return (
        (group["WIND"] >= ts_min_kt)
        & (group["WIND"] < hurricane_min_kt)
        & (group["distance_km"] < impact_radius_km - radius_reduction_km)
    )


def impact_mask(group, impact_radius_km, **kwargs):
    """Timesteps meeting the intensity + proximity criteria."""
    return (
        hurricane_mask(group, impact_radius_km,
                       kwargs.get("hurricane_min_kt", HURRICANE_MIN_KT))
        | ts_mask(group, impact_radius_km, **kwargs)
    )


def median_timestep_hours(group, default=6.0):
    """Median spacing between track points, in hours.

    IBTrACS v04 interleaves 3- and 6-hourly records, so duration metrics must
    not assume 6 hours.
    """
    if len(group) < 2:
        return default
    step = group["DATE"].diff().dt.total_seconds().median() / 3600
    return float(step) if np.isfinite(step) and step > 0 else default


# ── PER-STORM SUMMARY ─────────────────────────────────────────────────────────
def summarise_storm(sid, group, impact_radius_km, iso3, **criteria):
    """Summarise one storm's track into a row of metrics.

    ``group`` is all track points for ``sid`` (from :func:`load_ibtracs` plus
    :func:`add_distance`); ``**criteria`` is passed to the mask functions.
    Storms that never entered the buffer return a mostly-empty row.
    """
    group = group.sort_values("DATE").copy()
    dist, wind = group["distance_km"], group["WIND"]
    near = dist < impact_radius_km

    formation = group.iloc[0]
    formation_date = formation["DATE"]
    closest_row = group.loc[dist.idxmin()]
    closest_approach_date = closest_row["DATE"]

    if near.any():
        max_wind_kt = float(wind[near].max())
        dist_at_max = float(group.loc[wind[near].idxmax(), "distance_km"])
        first_hit = group.loc[near, "DATE"].min()
        travel_days = max((first_hit - formation_date).days, 0)
    else:
        max_wind_kt = dist_at_max = travel_days = np.nan

    # Intensifying on approach: wind strictly increasing over the 48 h before
    # closest approach (number of points derived from the actual timestep).
    step_h = median_timestep_hours(group)
    n_pre = max(int(round(48 / step_h)), 2)
    pre_window = group[group["DATE"] < closest_approach_date].tail(n_pre)
    intensifying = (
        bool(pre_window["WIND"].is_monotonic_increasing) if len(pre_window) >= 2 else None
    )

    h_mask = hurricane_mask(group, impact_radius_km,
                            criteria.get("hurricane_min_kt", HURRICANE_MIN_KT))
    t_mask = ts_mask(group, impact_radius_km, **criteria)

    season = formation["SEASON"]

    return pd.Series({
        "sid": sid,
        "iso3": iso3,
        "name": formation["NAME"],
        "year": int(season) if pd.notna(season) else np.nan,
        "formation_date": formation_date,
        "formation_month": formation_date.month,
        "formation_lat": formation["LAT"],
        "formation_lon": formation["LON"],
        "within_buffer": bool(near.any()),
        "formation_wind_kmh": knots_to_kmh(formation["WIND"]),
        "closest_approach_date": closest_approach_date,
        "closest_approach_lat": closest_row["LAT"],
        "closest_approach_lon": closest_row["LON"],
        "closest_dist_km": closest_row["distance_km"],
        "travel_days": travel_days,
        "max_wind_knots": max_wind_kt,
        "max_wind_kmh": knots_to_kmh(max_wind_kt),
        "ss_category": saffir_simpson_category(max_wind_kt, units="kt"),
        "dist_at_max_intensity_km": dist_at_max,
        "intensifying_on_approach": intensifying,
        "hours_as_hurricane": round(int(h_mask.sum()) * step_h, 1),
        "hours_as_ts": round(int(t_mask.sum()) * step_h, 1),
        "track_timestep_hours": step_h,
    })


def summarise_storms(track_points, impact_radius_km, iso3, **criteria):
    """Summarise every SID in ``track_points``."""
    rows = [
        summarise_storm(sid, grp, impact_radius_km, iso3, **criteria)
        for sid, grp in track_points.groupby("SID")
    ]
    return pd.DataFrame(rows).reset_index(drop=True)
