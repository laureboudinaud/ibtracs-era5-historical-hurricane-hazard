"""Shared plotting style and reusable figures for the hazard analysis.

Every ``plot_*`` function saves a PNG to ``output_dir`` and returns the figure.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# ── STYLE ─────────────────────────────────────────────────────────────────────
SS_COLORS = {
    "TS":    "#7b68ee",
    "Cat 1": "#FFE135",
    "Cat 2": "#FFA500",
    "Cat 3": "#FF6B1A",
    "Cat 4": "#E8390E",
    "Cat 5": "#CC0000",
    "Non-TC": "#9E9E9E",   # troughs / waves — grey and last, never on the SS scale
}
FALLBACK_COLOR = "#aaaaaa"
WINDOW_COLORS = {"3d": "#0C447C", "5d": "#378ADD", "7d": "#B5D4F4"}

PCT_VALS = [25, 50, 75]
PCT_STYLES = [":", "--", "-."]
PCT_COLORS = ["#888888", "#444444", "#111111"]


def set_style(font="Open Sans", context="talk", font_scale=0.75, style="whitegrid"):
    """Apply the shared seaborn style, falling back to the default font."""
    from matplotlib import font_manager

    sns.set_theme(style=style, context=context)
    sns.set_context(context, font_scale=font_scale)

    if font in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = font
    else:
        print(f"  [info] font '{font}' not installed — using matplotlib default")


# ── HELPERS ───────────────────────────────────────────────────────────────────
def normalise_cat(c):
    """Normalise a category label to 'TS' or 'Cat 1'…'Cat 5'.

    Accepts the string forms that survive a CSV round-trip ('3', '3.0',
    'Cat 3'). Anything unrecognised passes through, which is how 'Non-TC'
    survives.
    """
    c = str(c).strip()
    if c in ("Tropical Storm", "TS"):
        return "TS"
    if c.startswith("Cat "):
        return c
    try:
        return f"Cat {int(float(c))}"
    except ValueError:
        return c


def cat_color(cat):
    """Colour for a category label."""
    return SS_COLORS.get(normalise_cat(cat), FALLBACK_COLOR)


def save_fig(fig, name, output_dir="outputs"):
    """Save at 150 dpi with a tight bounding box (labels drawn outside the axes
    depend on it), creating ``output_dir`` if needed."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved → {path}")
    return path


def category_handles(cats=None):
    """Legend handles, always ordered TS → Cat 5 → Non-TC."""
    present = list(SS_COLORS) if cats is None else [c for c in SS_COLORS if c in set(cats)]
    return [mpatches.Patch(color=SS_COLORS[c], label=c) for c in present]


def percentile_handles():
    """Legend handles for the P25/P50/P75 reference lines."""
    return [mlines.Line2D([], [], color=pc, lw=1.0, ls=ls, label=f"P{p}")
            for p, ls, pc in zip(PCT_VALS, PCT_STYLES, PCT_COLORS)]


def _clean(arr):
    """Float array with NaNs dropped."""
    arr = np.asarray(arr, dtype=float)
    return arr[~np.isnan(arr)]


def _median_label(ax, y, text, color):
    """Median annotation just outside the right spine, in data-space y."""
    ax.text(1.005, y, text, transform=ax.get_yaxis_transform(),
            fontsize=7.5, color=color, ha="left", va="center",
            zorder=5, clip_on=False)


# ── EXCEEDANCE ────────────────────────────────────────────────────────────────
def plot_exceedance(r3, r5, r7, country_name, output_dir="outputs",
                    record_years=None, filename="7_exceedance_thresholds.png"):
    """Exceedance curves for the three windows, plus a 3-day threshold table.

    Uses the Weibull plotting position ``i/(n+1)``, so the largest observed
    event is not given a 0 % exceedance probability. With ``record_years`` the
    table reports a return period in years, otherwise a frequency in events.

    Pass only events drawn by a defined search (the IBTrACS set): a curated
    list of damaging events has no denominator.
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    ax = axes[0]
    for (lbl, color), arr in zip(WINDOW_COLORS.items(), [r3, r5, r7]):
        data = np.sort(_clean(arr))
        if data.size == 0:
            continue
        exceedance = 1 - np.arange(1, data.size + 1) / (data.size + 1)
        ax.plot(data, exceedance, color=color, lw=2, label=f"{lbl} total")
    ax.set_xlabel("Rainfall threshold (mm)")
    ax.set_ylabel("Probability of exceedance")
    ax.set_title("Exceedance Probability Curves", fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
    ax.legend()
    ax.yaxis.grid(True, alpha=0.4)
    ax.xaxis.grid(True, alpha=0.4)
    ax.spines[["top", "right"]].set_visible(False)

    ax2 = axes[1]
    data3d = _clean(r3)
    n = data3d.size
    rows = []
    for p in [50, 60, 70, 75, 80, 85, 90, 95]:
        t = np.percentile(data3d, p) if n else np.nan
        n_ex = int((data3d >= t).sum()) if n else 0
        if n_ex == 0:
            freq = "—"
        elif record_years:
            freq = f"~1 in {record_years / n_ex:.1f} yr"
        else:
            freq = f"~1 in {n / n_ex:.1f} events"
        rows.append([f"P{p}  ({t:.1f} mm)", f"{n_ex}/{n}",
                     f"{n_ex / n * 100:.0f}%" if n else "—", freq])

    ax2.axis("off")
    tbl = ax2.table(
        cellText=rows,
        colLabels=["3d threshold", "Events exceeding", "% of events",
                   "Return period" if record_years else "Frequency"],
        cellLoc="center", loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.8)
    for (r, _c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2c3e6b")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f0f4ff")
    ax2.set_title("3-day Rainfall Threshold Analysis", fontweight="bold", pad=15)

    plt.suptitle(f"Exceedance Analysis — {country_name} Storm Rainfall",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


# ── HEATMAP ───────────────────────────────────────────────────────────────────
def plot_heatmap(r3, r5, r7, gmax, peak, dates, cats, country_name,
                 output_dir="outputs", filename="8_heatmap_all_events.png"):
    """Normalised heatmap of every event's metrics, ranked by 3-day total.

    Over a single ERA5-Land cell ``gmax`` equals the 7-day mean and is not a
    second metric — check the pixel count before reading it as one.
    """
    hm_df = pd.DataFrame(
        np.column_stack([r3, r5, r7, gmax, peak]),
        columns=["3d total", "5d total", "7d total", "7d grid max", "peak hourly"],
    )

    span = (hm_df.max() - hm_df.min()).replace(0, np.nan)   # guard zero-range columns
    hm_norm = ((hm_df - hm_df.min()) / span).fillna(0.5)

    order = hm_df["3d total"].sort_values(ascending=False).index
    hm_norm = hm_norm.loc[order].reset_index(drop=True)
    hm_annot = hm_df.loc[order].round(1).reset_index(drop=True)
    row_labels = [f"{dates[i]}  {cats[i]}" for i in order]

    fig, ax = plt.subplots(figsize=(9, max(8, len(row_labels) * 0.28)))
    sns.heatmap(hm_norm, ax=ax, cmap="YlOrRd", linewidths=0.3, linecolor="white",
                annot=hm_annot, fmt=".1f", annot_kws={"size": 7},
                yticklabels=row_labels,
                cbar_kws={"label": "Normalised value (0–1)"})
    ax.set_title(f"All events ranked by 3d rainfall, annotated with actual values "
                 f"({country_name})", fontweight="bold", pad=12)
    ax.tick_params(axis="y", labelsize=7)
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


# ── TREND ─────────────────────────────────────────────────────────────────────
def plot_trend_over_time(years, values, colors, country_name,
                         ylabel="3-day rainfall (mm)", output_dir="outputs",
                         filename="9_trend_over_time.png", show_fit=True):
    """Event-level metric against year, with an optional OLS fit.

    A significant slope here is not evidence of a climate trend: the sample is
    conditioned on the proximity and intensity filters, and the reanalysis is
    not homogeneous across the pre-satellite era.
    """
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    if len(years) != len(values):
        raise ValueError(f"years ({len(years)}) and values ({len(values)}) must align")
    ok = ~np.isnan(years) & ~np.isnan(values)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(years, values, c=colors, s=60, zorder=3,
               edgecolors="white", linewidths=0.5)

    if show_fit and ok.sum() >= 3:
        m, b, r_val, p_val, _ = sp_stats.linregress(years[ok], values[ok])
        xs = np.linspace(years[ok].min(), years[ok].max(), 10)
        ax.plot(xs, m * xs + b, color="#333", lw=1.5, linestyle="--",
                label=f"Linear trend  (r={r_val:.2f}, p={p_val:.2f})")
        ax.legend(fontsize=9)

    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} over time — {country_name}", fontweight="bold")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


# ── PER-EVENT BARS ────────────────────────────────────────────────────────────
def plot_stacked_event_bars(r3, r5, r7, colors, xlabels, country_name,
                            cats=None, output_dir="outputs",
                            filename="2a_rainfall_bars.png", bar_width=0.45,
                            figsize=(20, 5.5)):
    """Per-event rainfall: 3-day base with 5- and 7-day increments stacked as
    lighter segments of the same category colour."""
    r3, r5, r7 = (np.asarray(a, dtype=float) for a in (r3, r5, r7))
    inc_5d, inc_7d = r5 - r3, r7 - r5
    x = np.arange(len(r3))

    fig, ax = plt.subplots(figsize=figsize)
    for i, col in enumerate(colors):
        bottom = 0.0
        for vals, alpha in [(r3, 1.0), (inc_5d, 0.6), (inc_7d, 0.3)]:
            ax.bar(x[i], vals[i], color=col, width=bar_width, alpha=alpha,
                   bottom=bottom, zorder=3)
            bottom += vals[i]

    for arr, color, lbl in [(r3, "#3D3D3D", "3d"), (r5, "#757575", "5d"),
                            (r7, "#BBBBBB", "7d")]:
        med = np.nanmedian(arr)
        ax.axhline(med, color=color, linewidth=0.8, linestyle="--", zorder=4)
        _median_label(ax, med, f" median {lbl}  {med:.1f}", color)

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("mm", fontsize=9)
    ax.set_title("Accumulated rainfall (3d + 5d + 7d increments)",
                 fontsize=10, pad=5, loc="left", color="dimgray")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    alpha_patches = [
        mpatches.Patch(color="gray", alpha=1.0, label="3-day total"),
        mpatches.Patch(color="gray", alpha=0.6, label="5-day increment"),
        mpatches.Patch(color="gray", alpha=0.3, label="7-day increment"),
    ]
    median_handle = [mlines.Line2D([], [], color="#3D3D3D", lw=0.8, ls="--",
                                   label="Median")]
    ax.legend(handles=category_handles(cats) + alpha_patches + median_handle,
              ncol=10, fontsize=8, frameon=False,
              loc="upper left", bbox_to_anchor=(0, 1.18))
    fig.suptitle(f"ERA5-Land rainfall — storm events near {country_name}",
                 fontsize=11, y=0.98)
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


def plot_event_bars(values, colors, xlabels, title, ylabel, country_name,
                    cats=None, output_dir="outputs", filename="event_bars.png",
                    bar_width=0.45, figsize=(20, 5.5)):
    """Single-metric bar chart per event, coloured by category, with a median line."""
    values = np.asarray(values, dtype=float)
    x = np.arange(len(values))

    fig, ax = plt.subplots(figsize=figsize)
    for i, col in enumerate(colors):
        ax.bar(x[i], values[i], color=col, width=bar_width, zorder=3)

    med = np.nanmedian(values)
    ax.axhline(med, color="#3D3D3D", linewidth=0.8, linestyle="--", zorder=4)
    _median_label(ax, med, f" median {med:.1f}", "#3D3D3D")

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10, pad=5, loc="left", color="dimgray")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    median_handle = [mlines.Line2D([], [], color="#3D3D3D", lw=0.8, ls="--",
                                   label="Median")]
    ax.legend(handles=category_handles(cats) + median_handle,
              ncol=10, fontsize=8, frameon=False,
              loc="upper left", bbox_to_anchor=(0, 1.18))
    fig.suptitle(f"ERA5-Land rainfall — storm events near {country_name}",
                 fontsize=11, y=0.98)
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


# ── DISTRIBUTIONS ─────────────────────────────────────────────────────────────
def plot_distributions(r3, r5, r7, country_name, output_dir="outputs",
                       filename="3_rainfall_distributions.png"):
    """Histogram + KDE + percentile lines per accumulation window."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (lbl, color), arr in zip(axes, WINDOW_COLORS.items(), [r3, r5, r7]):
        data = _clean(arr)
        if data.size == 0:
            ax.set_visible(False)
            continue
        ax.hist(data, bins=12, color=color, alpha=0.5, edgecolor="white", density=True)
        if data.size > 1 and np.ptp(data) > 0:
            # Clip the kernel at the data range: rainfall cannot be negative.
            kde_x = np.linspace(data.min(), data.max(), 200)
            ax.plot(kde_x, sp_stats.gaussian_kde(data)(kde_x), color=color, lw=2)
        for p, ls in zip(PCT_VALS, PCT_STYLES):
            val = np.percentile(data, p)
            ax.axvline(val, color="#333", ls=ls, lw=1.2, label=f"P{p}: {val:.1f} mm")
        ax.set_title(f"{lbl} accumulated rainfall", fontweight="bold")
        ax.set_xlabel("mm")
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    plt.suptitle(f"Rainfall Distributions — {country_name} Storm Events",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig


# ── RAINFALL VS CATEGORY / PREDICTOR ──────────────────────────────────────────
def plot_rainfall_vs_category(arr, cats, country_name, window_label="3d",
                              output_dir="outputs", filename=None, seed=0):
    """Jittered strip plot of rainfall by category, with medians and sample sizes.

    The point is that the two are only loosely related — 'Non-TC' events carry
    no wind hazard at all.
    """
    rng = np.random.default_rng(seed)          # reproducible jitter
    arr = np.asarray(arr, dtype=float)
    cats = np.asarray(cats)
    cat_order = [c for c in SS_COLORS if c in set(cats)]

    fig, ax = plt.subplots(figsize=(7, 4))
    for i, cat in enumerate(cat_order):
        vals = arr[cats == cat]
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   color=SS_COLORS.get(cat, FALLBACK_COLOR),
                   s=60, alpha=0.85, edgecolors="white", linewidths=0.5, zorder=3)
        if len(vals):
            ax.plot([i - 0.25, i + 0.25], [np.nanmedian(vals)] * 2,
                    color="#333", lw=1.5, zorder=4)

    ymin, ymax = ax.get_ylim()
    for i, cat in enumerate(cat_order):
        ax.text(i, ymin + (ymax - ymin) * 0.01, f"n={(cats == cat).sum()}",
                ha="center", fontsize=7, color="#666")

    ax.set_xticks(range(len(cat_order)))
    ax.set_xticklabels(cat_order)
    ax.set_title(f"{window_label} total rainfall by category — {country_name}\n"
                 "(rainfall hazard ≠ wind category)",
                 fontsize=10, pad=8, loc="left", color="dimgray")
    ax.set_ylabel("mm", fontsize=9)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    save_fig(fig, filename or f"4_rainfall_vs_category_{window_label}.png", output_dir)
    plt.show()
    return fig


def plot_rainfall_vs(x, arrays, cats, labels, xlabel, country_name,
                     annotations=None, output_dir="outputs",
                     filename="5_rainfall_vs.png", figsize=None):
    """Scatter rainfall against a continuous predictor, one panel per window.

    ``x`` is the predictor (distance, wind speed…), ``arrays`` one rainfall
    series per panel, ``annotations`` optional per-point text such as the year.
    """
    x = np.asarray(x, dtype=float)
    cats = np.asarray(cats)
    fig, axes = plt.subplots(1, len(arrays), figsize=figsize or (6 * len(arrays), 5))
    axes = np.atleast_1d(axes)

    for ax, label, arr in zip(axes, labels, arrays):
        arr = np.asarray(arr, dtype=float)
        for cat in dict.fromkeys(cats):
            m = cats == cat
            ax.scatter(x[m], arr[m], color=SS_COLORS.get(cat, FALLBACK_COLOR),
                       label=cat, s=70, alpha=0.85,
                       edgecolors="white", linewidths=0.5, zorder=3)
        if annotations is not None:
            for xi, yi, txt in zip(x, arr, annotations):
                ax.annotate(str(txt), xy=(xi, yi), fontsize=6, color="#444444",
                            xytext=(4, 4), textcoords="offset points")
        clean = _clean(arr)
        for p, ls, pc in zip(PCT_VALS, PCT_STYLES, PCT_COLORS):
            val = np.percentile(clean, p)
            ax.axhline(val, color=pc, linewidth=1.0, linestyle=ls, zorder=2)
            ax.text(np.nanmax(x), val, f" P{p} {val:.0f} mm",
                    fontsize=7, color=pc, va="bottom", ha="left", zorder=4)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(f"{label} rainfall (mm)", fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)
        ax.xaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.3, zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)

    handles, labels_leg = axes[0].get_legend_handles_labels()
    pct = percentile_handles()
    fig.legend(handles + pct, labels_leg + [h.get_label() for h in pct],
               ncol=len(handles) + len(pct), fontsize=8, frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(f"Rainfall vs {xlabel.lower()} — {country_name}", fontsize=12, y=1.01)
    plt.tight_layout()
    save_fig(fig, filename, output_dir)
    plt.show()
    return fig
