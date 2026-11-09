"""Paper figures from the high-speed camera records of selected plug and slug runs.

Reads the frames and space-time arrays extracted by scripts/extract_camera_frames.py
(data/bronze_layer/camera/, local only) and writes, for each selected run, one PDF per panel:

    fig_camera_frames_<tag>.pdf      a few frames stacked, flow from left to right
    fig_camera_spacetime_<tag>.pdf   height above the bottom wall vs camera time at the middle column

Camera time t_cam starts at the first recorded frame of each video. The camera and the confocal
sensor were not synchronized, so t_cam = 0 is not the confocal t = 0.

Writes a new folder data/gold_layer/camera_figures_<N>/ (never overwrites) with run_manifest.json.

Usage: python scripts/camera_figures.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.flow_pattern_criteria import CHANNEL_SIDE                       # noqa: E402
from src.flow_pattern_report import (COL_W, PAPER_RC, T_CAMERA_LABEL, run_file_tag,   # noqa: E402
                                     run_label)
from src.output_versions import next_version_dir, write_manifest        # noqa: E402

CAMERA = ROOT / "data/bronze_layer/camera"

# runs of the confocal profiles of the manuscript; frame times [ms] chosen by inspection of the contact sheets
SELECTION = {
    "0305_06": [64.0, 238.0],     # plug: bubble nose, bubble tail (film below the bubble in both)
    "0314_06": [104.0, 138.0],    # slug: aerated liquid slug, bubble tail with the aerated region behind it
}


def walls(folder):
    """Rows of the inner top and bottom walls (dark, static lines of the mean image) and pixels per mm."""
    z = np.load(folder / "spacetime.npz")
    prof, sprof = np.median(z["mean_img"], axis=1), np.median(z["std_img"], axis=1)
    h = len(prof)
    static = sprof < np.percentile(sprof, 10) * 1.6
    p = np.where(static, prof, np.inf)
    top, bot = int(np.argmin(p[: h // 3])), 2 * h // 3 + int(np.argmin(p[2 * h // 3:]))
    return top, bot, (bot - top) / CHANNEL_SIDE


def frames_figure(folder, meta, times, out):
    files = sorted((folder / "frames").glob("f*_t*ms.png"))
    t_files = np.array([float(f.stem.split("_t")[1][:-2]) for f in files])
    h, w = meta["crop_shape"]
    fig, axs = plt.subplots(len(times), 1, figsize=(COL_W, len(times) * (COL_W * h / w + 0.2) + 0.35),
                            squeeze=False, layout="constrained")
    used = []
    for ax, tq in zip(axs[:, 0], times):
        k = int(np.argmin(np.abs(t_files - tq)))
        used.append(float(t_files[k]))
        ax.imshow(np.asarray(Image.open(files[k])), cmap="gray", vmin=0, vmax=255, aspect="equal")
        ax.set_axis_off()
        ax.set_title(f"$t_{{cam}}$ = {t_files[k]:.0f} ms", fontsize=7, pad=1.5, loc="left")
    axs[0, 0].annotate("", xy=(0.99, 1.06), xytext=(0.84, 1.06), xycoords="axes fraction",
                       arrowprops=dict(arrowstyle="-|>", color="0.2", lw=0.8))
    axs[0, 0].text(0.83, 1.06, "Flow", transform=axs[0, 0].transAxes, ha="right", va="center", fontsize=7)
    fig.suptitle(run_label(meta["pattern"], meta["jl_m_s"], meta["jg_m_s"]), fontsize=7)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return used


def spacetime_figure(folder, meta, out):
    z = np.load(folder / "spacetime.npz")
    top, bot, px_mm = walls(folder)
    mid = len(z["probe_cols"]) // 2
    img = z["probes"][mid].T.astype(np.float32)          # rows x frames
    d = meta["display"]
    x = np.clip((img - d["lo"]) / (d["hi"] - d["lo"]), 0, 1) ** (1 / d["gamma"])
    t_ms = z["t_s"] * 1e3
    rows = img.shape[0]
    fig, ax = plt.subplots(figsize=(COL_W, 1.45), layout="constrained")
    ax.imshow(x, cmap="gray", aspect="auto", vmin=0, vmax=1,
              extent=[t_ms[0], t_ms[-1], (bot - rows) / px_mm, bot / px_mm], origin="upper")
    ax.set_ylim(-0.3, CHANNEL_SIDE + 0.3)
    ax.set_xlabel(T_CAMERA_LABEL)
    ax.set_ylabel("Height [mm]")
    ax.grid(False)
    ax.set_title(run_label(meta["pattern"], meta["jl_m_s"], meta["jg_m_s"]), fontsize=7, loc="left")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return {"top_row": top, "bottom_row": bot, "px_per_mm": px_mm}


def main():
    plt.rcParams.update(PAPER_RC)
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False})
    out = next_version_dir(ROOT / "data/gold_layer", "camera_figures")
    record = {}
    for di, times in SELECTION.items():
        folder = next(CAMERA.glob(f"run_*_{di}"))
        meta = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
        tag = run_file_tag(meta["pattern"], meta["jl_m_s"], meta["jg_m_s"])
        used = frames_figure(folder, meta, times, out / f"fig_camera_frames_{tag}.pdf")
        scale = spacetime_figure(folder, meta, out / f"fig_camera_spacetime_{tag}.pdf")
        record[tag] = {"run_id": meta["run_id"], "date_identifier": di, "fps": meta["fps"],
                       "frame_times_ms": used, **scale}
        print(tag, used, scale)
    write_manifest(out, "Camera frames and space-time images of selected plug and slug runs (camera time, not "
                        "synchronized with the confocal record)", runs=record)


if __name__ == "__main__":
    main()
