"""Slug-flow quantities of the plug and slug runs (Fluids paper, Stage 7 exploration).

Confocal (silver layer, processed_v2): slug units (src/slug_units.py) -> slug frequency, slug and
bubble durations, intermittency, length ratio, lengths with U = j_m, film below the bubbles.
Camera (bronze layer, scripts/extract_camera_frames.py): preliminary velocity of the interfaces.

Writes a new folder data/gold_layer/slug_analysis_<N>/ (never overwrites) with CSV tables,
run_manifest.json and the figures: PDF only, one panel per file, runs identified in the figures by
pattern and superficial velocities (run ids stay in the CSV tables for traceability).

Usage: python scripts/slug_flow_analysis.py
"""
import itertools
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.flow_pattern_criteria import CHANNEL_SIDE, load_run          # noqa: E402
from src.flow_pattern_report import (COL_W, PAPER_RC, SCATTER_H, T_CONFOCAL_LABEL,   # noqa: E402
                                     run_file_tag, run_label)
from src.output_versions import next_version_dir, write_manifest     # noqa: E402
from src.slug_units import (SLUG_PARAMS, SLUG_PARAM_GRID, detect_slug_units,   # noqa: E402
                            film_below_bubbles, unit_statistics)

CRITERIA_VERSION = "objective_criteria_11"
RUNS_TABLE = ROOT / "data/gold_layer" / CRITERIA_VERSION / "processed_v2/objective_flow_pattern.pkl"
SILVER_V2 = ROOT / "data/silver_layer/confocal_runs_processed_v2"
SILVER_RAW = ROOT / "data/silver_layer/confocal_runs"
CAMERA = ROOT / "data/bronze_layer/camera"

COLORS = {"plug": "#2a78d6", "slug": "#eb6834"}     # reference categorical slots 1 and 2
MARKERS = {"plug": "o", "slug": "^"}
FILM_SERIES = [34, 33, 32, 31, 30, 28]               # plug, j_g ~ 0.52 m/s, j_l 0.12 -> 0.79 m/s (kHz runs)
QC_RUNS = [33, 28, 56, 67, 25, 112]
CAMERA_SPREAD_MAX = 0.10                             # max relative spread of the three probe-pair velocities
WINDOW_S = 5.0                                       # window for the steadiness of the slug frequency [s]
WINDOW_MIN_RECORD = 20.0                             # records shorter than this give no window statistics [s]

plt.rcParams.update(PAPER_RC)
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "0.9"})


def save(fig, folder, name):
    """Figures are written as PDF only."""
    fig.savefig(folder / f"{name}.pdf")
    plt.close(fig)


def scatter(ax, df, x, y, flag=None, **kw):
    """Plug and slug markers; runs with df[flag] True are drawn hollow (uncertain)."""
    for pat, g in df.groupby("pattern"):
        bad = g[flag].astype(bool) if flag else pd.Series(False, index=g.index)
        ax.scatter(g.loc[~bad, x], g.loc[~bad, y], s=22, marker=MARKERS[pat], color=COLORS[pat],
                   edgecolor="white", linewidth=0.6, label=pat, zorder=3, **kw)
        if bad.any():
            ax.scatter(g.loc[bad, x], g.loc[bad, y], s=22, marker=MARKERS[pat], facecolor="none",
                       edgecolor=COLORS[pat], linewidth=0.8, label=f"{pat} (uncertain)", zorder=3, **kw)


def window_frequency(units, duration):
    """Steadiness of the slug frequency over the record: slugs starting in consecutive WINDOW_S windows.

    Returns the minimum, maximum and coefficient of variation of the window frequencies (NaN for records
    shorter than WINDOW_MIN_RECORD). Counts the slugs of the complete units (the last slug is not included).
    """
    if duration < WINDOW_MIN_RECORD or units.empty:
        return dict.fromkeys(["f_win_min", "f_win_max", "f_win_cv"], np.nan)
    starts = np.arange(0.0, duration - WINDOW_S + 1e-9, WINDOW_S)
    f = np.array([((units.s_start >= a) & (units.s_start < a + WINDOW_S)).sum() / WINDOW_S for a in starts])
    return {"f_win_min": f.min(), "f_win_max": f.max(), "f_win_cv": f.std() / f.mean() if f.mean() > 0 else np.nan}


# ---------------------------------------------------------------- camera (preliminary)
def camera_velocity(folder):
    """Interface velocity from the upper-half signal of the three probe columns (cross-correlation)."""
    z = np.load(folder / "spacetime.npz")
    prof, sprof = np.median(z["mean_img"], axis=1), np.median(z["std_img"], axis=1)
    h = len(prof)
    static = sprof < np.percentile(sprof, 10) * 1.6          # walls: dark and static
    p = np.where(static, prof, np.inf)
    top, bot = int(np.argmin(p[: h // 3])), 2 * h // 3 + int(np.argmin(p[2 * h // 3:]))
    px_mm = (bot - top) / CHANNEL_SIDE
    t, cols = z["t_s"], z["probe_cols"]
    dt = float(np.median(np.diff(t)))
    sig = z["probes"].astype(np.float32)[:, :, top + 10:(top + bot) // 2].mean(axis=2)

    def xvel(a, b, dx):
        a, b = np.diff(a), np.diff(b)
        a, b = (a - a.mean()) / a.std(), (b - b.mean()) / b.std()
        n = len(a)
        lags = np.arange(max(1, int(dx / 10 / dt)), min(int(dx / 0.2 / dt), n // 2) + 1)
        cc = np.array([np.mean(a[: n - L] * b[L:]) for L in lags])
        k, off = int(np.argmax(cc)), 0.0
        if 0 < k < len(cc) - 1:
            y0, y1, y2 = cc[k - 1: k + 2]
            den = y0 - 2 * y1 + y2
            off = 0.5 * (y0 - y2) / den if den else 0.0
        return dx / ((lags[k] + off) * dt), float(cc[k])

    pairs = [(0, 2), (0, 1), (1, 2)]
    res = [xvel(sig[i], sig[j], (cols[j] - cols[i]) / px_mm * 1e-3) for i, j in pairs]
    u = [r[0] for r in res]
    return {"px_mm": px_mm, "U_cam": u[0], "cc_cam": res[0][1], "U_spread": (max(u) - min(u)) / u[0]}


def main():
    warnings.filterwarnings("ignore")
    runs = pd.read_pickle(RUNS_TABLE)
    runs = runs[runs.reduced_pattern.isin(["plug", "slug"])].copy()
    runs["pattern"] = runs.reduced_pattern

    out = next_version_dir(ROOT / "data/gold_layer", "slug_analysis")
    print("output:", out)

    grid = [dict(zip(SLUG_PARAM_GRID, v)) for v in itertools.product(*SLUG_PARAM_GRID.values())]
    rows, units_all, sens, data = [], [], [], {}
    for r in runs.itertuples():
        df = load_run(int(r.run_id), SILVER_V2, raw_path=SILVER_RAW)
        data[r.run_id] = df
        jm = r.jl + r.jg
        units, info = detect_slug_units(df, jm)
        st = unit_statistics(units, info, r.jl, r.jg)
        f_grid = []
        for g in grid:
            u2, i2 = detect_slug_units(df, jm, g)
            f_grid.append(i2["n_slugs"] / i2["duration"])
            sens.append({"run_id": r.run_id, **g, "f_slug": f_grid[-1]})
        rows.append({"run_id": r.run_id, "date_identifier": r.date_identifier, "pattern": r.pattern,
                     "jl": r.jl, "jg": r.jg, "sampling_hz": r.sampling_hz, "f_undetected": r.f_undetected,
                     "f_bridge": r.bridge_freq_hz, **st, **window_frequency(units, info["duration"]),
                     "dt": info["dt"], "t_bubble_min": info["t_bubble_min"],
                     "f_slug_p10": np.percentile(f_grid, 10), "f_slug_p90": np.percentile(f_grid, 90)})
        units_all.append(units.assign(run_id=r.run_id))
    res = pd.DataFrame(rows)
    res["lambda_l"] = res.jl / res.jm
    res["rate"] = np.where(res.sampling_hz > 200, "kHz", "0.1 kHz")
    res["f_uncertain"] = res.f_slug_p90 > 2 * res.f_slug_p10          # frequency depends on the parameters
    res["ts_resolution"] = res.t_s_med <= 1.5 * res.dt                 # median slug of one sample
    res["ratio_uncertain"] = res.f_uncertain | res.ts_resolution        # t_b / t_s depends on both

    cam = []
    for folder in sorted(CAMERA.glob("run_*")):
        meta = json.loads((folder / "metadata.json").read_text())
        cam.append({"run_id": meta["run_id"], **camera_velocity(folder)})
    res = res.merge(pd.DataFrame(cam), on="run_id", how="left")
    res["U_cam_ok"] = res.U_spread < CAMERA_SPREAD_MAX

    res.to_csv(out / "slug_runs.csv", index=False)
    pd.concat(units_all).to_csv(out / "slug_units.csv", index=False)
    pd.DataFrame(sens).to_csv(out / "slug_sensitivity.csv", index=False)

    # ---------------------------------------------------------------- figures (PDF only, one panel each)
    one = (COL_W, SCATTER_H)
    res["jg_jm"] = res.jg / res.jm
    res["jg_jl"] = res.jg / res.jl

    # slug frequency vs j_l and vs j_g, with the sensitivity band
    for x, lab in [("jl", "$j_l$ [m/s]"), ("jg", "$j_g$ [m/s]")]:
        fig, ax = plt.subplots(figsize=one, layout="constrained")
        for pat, g in res.groupby("pattern"):
            ax.errorbar(g[x], g.f_slug, yerr=[g.f_slug - g.f_slug_p10, g.f_slug_p90 - g.f_slug], fmt="none",
                        ecolor=COLORS[pat], elinewidth=0.8, alpha=0.6, zorder=2)
        scatter(ax, res, x, "f_slug", flag="f_uncertain")
        ax.set_xlabel(lab)
        ax.set_ylabel("Slug frequency $f_s$ [Hz]")
        ax.legend(frameon=False)
        save(fig, out, f"fig_slug_frequency_vs_{x}")

    # slug frequency vs liquid-bridge frequency
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    scatter(ax, res, "f_bridge", "f_slug", flag="f_uncertain")
    lim = [0, max(res.f_slug.max(), res.f_bridge.max()) * 1.05]
    ax.plot(lim, lim, color="0.4", lw=0.8, ls="--", zorder=1)
    ax.text(lim[1] * 0.62, lim[1] * 0.55, "1:1", color="0.4")
    ax.set_xlabel("Liquid bridge frequency $f_{bridge}$ [Hz]")
    ax.set_ylabel("Slug frequency $f_s$ [Hz]")
    ax.legend(frameon=False)
    save(fig, out, "fig_slug_vs_bridge_frequency")

    # intermittency vs gas fraction of the mixture velocity (unit at j_m, slug without gas: j_g / j_m)
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    scatter(ax, res, "jg_jm", "intermittency", flag="f_uncertain")
    ax.plot([0, 1], [0, 1], color="0.4", lw=0.8, ls="--", zorder=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$j_g/j_m$ [-]")
    ax.set_ylabel(r"Intermittency $\gamma$ [-]")
    ax.legend(frameon=False, loc="upper left")
    save(fig, out, "fig_intermittency")

    # lengths with U = j_m and the measured ratio
    for y, flag, lab, name in [("Ls_Dh", "ts_resolution", "$L_S/D_h$ ($U = j_m$) [-]", "fig_slug_length"),
                               ("Lb_Dh", "f_uncertain", "$L_B/D_h$ ($U = j_m$) [-]", "fig_bubble_length")]:
        fig, ax = plt.subplots(figsize=one, layout="constrained")
        scatter(ax, res, "jm", y, flag=flag)
        ax.set_xlabel("$j_m$ [m/s]")
        ax.set_ylabel(lab)
        ax.set_yscale("log")
        ax.legend(frameon=False)
        save(fig, out, name)
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    scatter(ax, res, "jg_jl", "ratio_bs", flag="ratio_uncertain")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$j_g/j_l$ [-]")
    ax.set_ylabel("Median $t_b/t_s$ [-]")
    ax.legend(frameon=False)
    save(fig, out, "fig_length_ratio")

    # film below the elongated bubbles, plug series at j_g ~ 0.52 m/s
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    blues = plt.get_cmap("Blues")(np.linspace(0.35, 0.95, len(FILM_SERIES)))
    film_rows = []
    for c, rid in zip(blues, FILM_SERIES):
        units = pd.concat(units_all).query("run_id == @rid")
        x, film, nb = film_below_bubbles(data[rid], units)
        rr = res.set_index("run_id").loc[rid]
        ax.plot(x, film / CHANNEL_SIDE, color=c, lw=1.6, label=f"$j_l$ = {rr.jl:.2f} m/s")
        film_rows += [{"run_id": rid, "x": xi, "film_mm": fi, "n_bubbles": ni} for xi, fi, ni in zip(x, film, nb)]
    ax.set_xlabel("Position in the bubble $x_b$ [-]")
    ax.set_ylabel(r"Bottom film $\delta_f/D_h$ [-]")
    ax.legend(frameon=False, title="plug, $j_g$ = 0.52 m/s", fontsize=7)
    save(fig, out, "fig_film_below_bubbles_profiles")
    fr = pd.DataFrame(film_rows)
    fr.to_csv(out / "film_below_bubbles.csv", index=False)

    agg = fr.groupby("run_id").apply(lambda g: pd.Series({
        "front": g.film_mm[g.x < 0.15].median(), "end": g.film_mm[g.x > 0.85].median()})).reset_index()
    agg = agg.merge(res[["run_id", "jl"]], on="run_id")
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    ax.plot(agg.jl, agg.front / CHANNEL_SIDE, "o-", color=COLORS["plug"], label="Front of the bubble ($x_b < 0.15$)")
    ax.plot(agg.jl, agg.end / CHANNEL_SIDE, "s--", color="#1c5cab", label="End of the bubble ($x_b > 0.85$)")
    ax.set_xlabel("$j_l$ [m/s]")
    ax.set_ylabel(r"Bottom film $\delta_f/D_h$ [-]")
    ax.legend(frameon=False, fontsize=7, title="plug, $j_g$ = 0.52 m/s")
    save(fig, out, "fig_film_below_bubbles_front_end")

    # camera interface velocity vs j_m (preliminary)
    fig, ax = plt.subplots(figsize=one, layout="constrained")
    ok = res[res.U_cam_ok]
    scatter(ax, ok, "jm", "U_cam")
    bad = res[~res.U_cam_ok & res.U_cam.notna()]
    ax.scatter(bad.jm, bad.U_cam.clip(upper=8), s=12, facecolor="none", edgecolor="0.6", lw=0.6,
               label="Probe pairs disagree", zorder=2)
    lim = [0, max(4, ok.jm.max() * 1.1)]
    ax.plot(lim, lim, color="0.4", lw=0.8, ls="--", zorder=1)
    ax.text(lim[1] * 0.75, lim[1] * 0.68, "$U = j_m$", color="0.4")
    ax.set_xlim(lim)
    ax.set_ylim(0, 8)
    ax.set_xlabel("$j_m$ [m/s]")
    ax.set_ylabel("Camera interface velocity $U$ [m/s]")
    ax.legend(frameon=False, fontsize=7)
    save(fig, out, "fig_camera_velocity")

    # quality control: detection on the confocal traces, one figure per run
    qc = out / "detection_qc"
    qc.mkdir()
    for rid in QC_RUNS:
        df = data[rid]
        rr = res.set_index("run_id").loc[rid]
        t0 = 0.0 if rr.sampling_hz > 200 else 5.0
        t1 = t0 + (1.0 if rr.sampling_hz > 200 else 10.0)
        w = (df.time >= t0) & (df.time < t1)
        mw = df.mixer_wave.astype(bool)
        fig, ax = plt.subplots(figsize=(COL_W, SCATTER_H), layout="constrained")
        ax.plot(df.time[w & ~mw], df.t2[w & ~mw], ".", ms=1.5, color="0.15", label=r"$\delta_2$")
        ax.plot(df.time[w & mw], np.full((w & mw).sum(), 8.3), "|", ms=4, color=COLORS["slug"],
                label="Removed as mixer wave")
        units = pd.concat(units_all).query("run_id == @rid")
        shown = False
        for u in units.itertuples():
            if u.s_end > t0 and u.s_start < t1:
                ax.axvspan(u.s_start, u.s_end, color=COLORS["plug"], alpha=0.18, lw=0,
                           label=None if shown else "Detected liquid slug")
                shown = True
        ax.set_xlim(t0, t1)
        ax.set_ylim(-0.2, 9.4)
        ax.set_ylabel(r"$\delta_2$ [mm]")
        ax.set_xlabel(T_CONFOCAL_LABEL)
        ax.set_title(f"{run_label(rr.pattern, rr.jl, rr.jg)}, {rr.sampling_hz / 1000:.1f} kHz", fontsize=7, loc="left")
        ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncols=3,
                  handlelength=1.2, columnspacing=1.0)
        save(fig, qc, f"fig_detection_qc_{run_file_tag(rr.pattern, rr.jl, rr.jg)}")

    write_manifest(out, "Slug units from the confocal signal (plug and slug runs) and preliminary camera "
                        "interface velocity", criteria_version=CRITERIA_VERSION, slug_params=SLUG_PARAMS,
                   slug_param_grid=SLUG_PARAM_GRID, camera_spread_max=CAMERA_SPREAD_MAX,
                   film_series=FILM_SERIES, window_s=WINDOW_S, window_min_record=WINDOW_MIN_RECORD,
                   n_runs=len(res))
    print(res.groupby("pattern")[["f_slug", "f_bridge", "intermittency", "ratio_bs", "Ls_Dh", "Lb_Dh"]].median().round(3))
    print("camera velocity ok:", int(res.U_cam_ok.sum()), "of", int(res.U_cam.notna().sum()))


if __name__ == "__main__":
    main()
