"""Slug-flow quantities of the plug and slug runs (Fluids paper, Stage 7 exploration).

Confocal (silver layer, processed_v2): slug units (src/slug_units.py) -> slug frequency, slug and
bubble durations, intermittency, length ratio, lengths with U = j_m, film below the bubbles.
Camera (bronze layer, scripts/extract_camera_frames.py): preliminary velocity of the interfaces.

Writes a new folder data/gold_layer/slug_analysis_<N>/ (never overwrites) with CSV tables,
figures and run_manifest.json.

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

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "0.9", "grid.linewidth": 0.6})


def scatter(ax, df, x, y, flag=None, **kw):
    """Plug and slug markers; runs with df[flag] True are drawn hollow (uncertain)."""
    for pat, g in df.groupby("pattern"):
        bad = g[flag].astype(bool) if flag else pd.Series(False, index=g.index)
        ax.scatter(g.loc[~bad, x], g.loc[~bad, y], s=22, marker=MARKERS[pat], color=COLORS[pat],
                   edgecolor="white", linewidth=0.6, label=pat, zorder=3, **kw)
        if bad.any():
            ax.scatter(g.loc[bad, x], g.loc[bad, y], s=22, marker=MARKERS[pat], facecolor="none",
                       edgecolor=COLORS[pat], linewidth=0.8, label=f"{pat} (uncertain)", zorder=3, **kw)


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
                     "f_bridge": r.bridge_freq_hz, **st,
                     "dt": info["dt"], "t_bubble_min": info["t_bubble_min"],
                     "f_slug_p10": np.percentile(f_grid, 10), "f_slug_p90": np.percentile(f_grid, 90)})
        units_all.append(units.assign(run_id=r.run_id))
    res = pd.DataFrame(rows)
    res["lambda_l"] = res.jl / res.jm
    res["rate"] = np.where(res.sampling_hz > 200, "kHz", "0.1 kHz")
    res["f_uncertain"] = res.f_slug_p90 > 2 * res.f_slug_p10          # frequency depends on the parameters
    res["ts_resolution"] = res.t_s_med <= 1.5 * res.dt                 # median slug of one sample

    cam = []
    for folder in sorted(CAMERA.glob("run_*")):
        meta = json.loads((folder / "metadata.json").read_text())
        cam.append({"run_id": meta["run_id"], **camera_velocity(folder)})
    res = res.merge(pd.DataFrame(cam), on="run_id", how="left")
    res["U_cam_ok"] = res.U_spread < CAMERA_SPREAD_MAX

    res.to_csv(out / "slug_runs.csv", index=False)
    pd.concat(units_all).to_csv(out / "slug_units.csv", index=False)
    pd.DataFrame(sens).to_csv(out / "slug_sensitivity.csv", index=False)

    # ---------------------------------------------------------------- figures
    # 1 slug frequency vs j_l and j_g, with the sensitivity band
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True, layout="constrained")
    for ax, x, lab in zip(axs, ["jl", "jg"], ["$j_l$ [m/s]", "$j_g$ [m/s]"]):
        for pat, g in res.groupby("pattern"):
            ax.errorbar(g[x], g.f_slug, yerr=[g.f_slug - g.f_slug_p10, g.f_slug_p90 - g.f_slug], fmt="none",
                        ecolor=COLORS[pat], elinewidth=0.8, alpha=0.6, zorder=2)
        scatter(ax, res, x, "f_slug", flag="f_uncertain")
        ax.set_xlabel(lab)
    axs[0].set_ylabel("slug frequency $f_s$ [Hz]")
    axs[0].legend(frameon=False)
    fig.savefig(out / "fig_slug_frequency.png", dpi=200)
    fig.savefig(out / "fig_slug_frequency.pdf")
    plt.close(fig)

    # 2 slug frequency vs liquid-bridge frequency
    fig, ax = plt.subplots(figsize=(3.4, 3.0), layout="constrained")
    scatter(ax, res, "f_bridge", "f_slug", flag="f_uncertain")
    lim = [0, max(res.f_slug.max(), res.f_bridge.max()) * 1.05]
    ax.plot(lim, lim, color="0.4", lw=0.8, ls="--", zorder=1)
    ax.text(lim[1] * 0.62, lim[1] * 0.55, "1:1", color="0.4")
    ax.set_xlabel("liquid-bridge frequency $f_{bridge}$ [Hz]")
    ax.set_ylabel("slug frequency $f_s$ [Hz]")
    ax.legend(frameon=False)
    fig.savefig(out / "fig_slug_vs_bridge_frequency.png", dpi=200)
    plt.close(fig)

    # 3 intermittency vs gas fraction of the mixture velocity (no-slip, liquid slug without gas: beta = j_g / j_m)
    fig, ax = plt.subplots(figsize=(3.4, 3.0), layout="constrained")
    res["jg_jm"] = res.jg / res.jm
    scatter(ax, res, "jg_jm", "intermittency", flag="f_uncertain")
    ax.plot([0, 1], [0, 1], color="0.4", lw=0.8, ls="--", zorder=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$j_g / j_m$ [-]")
    ax.set_ylabel(r"intermittency $\sum t_b / \sum(t_s + t_b)$ [-]")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(out / "fig_intermittency.png", dpi=200)
    plt.close(fig)

    # 4 lengths with U = j_m and the measured ratio
    fig, axs = plt.subplots(1, 3, figsize=(8.5, 2.8), layout="constrained")
    scatter(axs[0], res, "jm", "Ls_Dh", flag="ts_resolution")
    axs[0].set_ylabel("$L_S / D_h$ (U = $j_m$) [-]")
    scatter(axs[1], res, "jm", "Lb_Dh", flag="f_uncertain")
    axs[1].set_ylabel("$L_B / D_h$ (U = $j_m$) [-]")
    for ax in axs[:2]:
        ax.set_xlabel("$j_m$ [m/s]")
        ax.set_yscale("log")
    res["jg_jl"] = res.jg / res.jl
    scatter(axs[2], res, "jg_jl", "ratio_bs", flag="ts_resolution")
    axs[2].set_xscale("log")
    axs[2].set_yscale("log")
    axs[2].set_xlabel("$j_g / j_l$ [-]")
    axs[2].set_ylabel("median $t_b / t_s$ = $L_B / L_S$ [-]")
    axs[0].legend(frameon=False)
    fig.savefig(out / "fig_slug_lengths.png", dpi=200)
    plt.close(fig)

    # 5 film below the elongated bubbles, plug series at j_g ~ 0.52 m/s
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), layout="constrained")
    blues = plt.get_cmap("Blues")(np.linspace(0.35, 0.95, len(FILM_SERIES)))
    film_rows = []
    for c, rid in zip(blues, FILM_SERIES):
        units = pd.concat(units_all).query("run_id == @rid")
        x, film, nb = film_below_bubbles(data[rid], units)
        rr = res.set_index("run_id").loc[rid]
        axs[0].plot(x, film / CHANNEL_SIDE, color=c, lw=1.6, label=f"$j_l$ = {rr.jl:.2f}")
        film_rows += [{"run_id": rid, "x": xi, "film_mm": fi, "n_bubbles": ni} for xi, fi, ni in zip(x, film, nb)]
    axs[0].set_xlabel("position in the bubble $(t - t_{slug,end}) / t_b$ [-]")
    axs[0].set_ylabel(r"bottom film $\delta_f / D_h$ [-]")
    axs[0].legend(frameon=False, title="$j_g$ = 0.52 m/s", fontsize=7)
    fr = pd.DataFrame(film_rows)
    agg = fr.groupby("run_id").apply(lambda g: pd.Series({
        "front": g.film_mm[g.x < 0.15].median(), "end": g.film_mm[g.x > 0.85].median()})).reset_index()
    agg = agg.merge(res[["run_id", "jl"]], on="run_id")
    axs[1].plot(agg.jl, agg.front / CHANNEL_SIDE, "o-", color=COLORS["plug"], label="front of the bubble ($x$ < 0.15)")
    axs[1].plot(agg.jl, agg.end / CHANNEL_SIDE, "s--", color="#1c5cab", label="end of the bubble ($x$ > 0.85)")
    axs[1].set_xlabel("$j_l$ [m/s]")
    axs[1].set_ylabel(r"$\delta_f / D_h$ [-]")
    axs[1].legend(frameon=False, fontsize=7)
    fig.savefig(out / "fig_film_below_bubbles.png", dpi=200)
    plt.close(fig)
    fr.to_csv(out / "film_below_bubbles.csv", index=False)

    # 6 camera interface velocity vs j_m (preliminary)
    fig, ax = plt.subplots(figsize=(3.4, 3.0), layout="constrained")
    ok = res[res.U_cam_ok]
    scatter(ax, ok, "jm", "U_cam")
    bad = res[~res.U_cam_ok & res.U_cam.notna()]
    ax.scatter(bad.jm, bad.U_cam.clip(upper=8), s=12, facecolor="none", edgecolor="0.6", lw=0.6,
               label="probe pairs disagree", zorder=2)
    lim = [0, max(4, ok.jm.max() * 1.1)]
    ax.plot(lim, lim, color="0.4", lw=0.8, ls="--", zorder=1)
    ax.text(lim[1] * 0.75, lim[1] * 0.68, "$U = j_m$", color="0.4")
    ax.set_xlim(lim)
    ax.set_ylim(0, 8)
    ax.set_xlabel("$j_m$ [m/s]")
    ax.set_ylabel("camera interface velocity $U$ [m/s]")
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(out / "fig_camera_velocity.png", dpi=200)
    plt.close(fig)

    # 7 quality control: detection on the confocal traces
    fig, axs = plt.subplots(len(QC_RUNS), 1, figsize=(8.0, 1.5 * len(QC_RUNS)), layout="constrained")
    for ax, rid in zip(axs, QC_RUNS):
        df = data[rid]
        rr = res.set_index("run_id").loc[rid]
        t0 = 0.0 if rr.sampling_hz > 200 else 5.0
        t1 = t0 + (1.0 if rr.sampling_hz > 200 else 10.0)
        w = (df.time >= t0) & (df.time < t1)
        mw = df.mixer_wave.astype(bool)
        ax.plot(df.time[w & ~mw], df.t2[w & ~mw], ".", ms=1.5, color="0.15")
        ax.plot(df.time[w & mw], np.full((w & mw).sum(), 8.3), "|", ms=4, color="#eb6834", label="flagged as mixer wave")
        units = pd.concat(units_all).query("run_id == @rid")
        for u in units.itertuples():
            if u.s_end > t0 and u.s_start < t1:
                ax.axvspan(u.s_start, u.s_end, color="#2a78d6", alpha=0.18, lw=0)
        ax.set_xlim(t0, t1)
        ax.set_ylim(-0.2, 8.6)
        ax.set_ylabel("$\\delta_2$ [mm]")
        ax.set_title(f"run {rid} ({rr.date_identifier}), {rr.pattern}, $j_l$ = {rr.jl:.2f}, $j_g$ = {rr.jg:.2f} m/s, "
                     f"{rr.sampling_hz / 1000:.2f} kHz: $f_s$ = {rr.f_slug:.2f} Hz ({rr.f_slug_p10:.2f}-{rr.f_slug_p90:.2f}), "
                     f"$f_{{bridge}}$ = {rr.f_bridge:.2f} Hz", fontsize=7, loc="left")
    axs[-1].set_xlabel("t [s]  (shaded: detected liquid slugs)")
    fig.savefig(out / "fig_detection_qc.png", dpi=170)
    plt.close(fig)

    write_manifest(out, "Slug units from the confocal signal (plug and slug runs) and preliminary camera "
                        "interface velocity", criteria_version=CRITERIA_VERSION, slug_params=SLUG_PARAMS,
                   slug_param_grid=SLUG_PARAM_GRID, camera_spread_max=CAMERA_SPREAD_MAX,
                   film_series=FILM_SERIES, n_runs=len(res))
    print(res.groupby("pattern")[["f_slug", "f_bridge", "intermittency", "ratio_bs", "Ls_Dh", "Lb_Dh"]].median().round(3))
    print("camera velocity ok:", int(res.U_cam_ok.sum()), "of", int(res.U_cam.notna().sum()))


if __name__ == "__main__":
    main()
