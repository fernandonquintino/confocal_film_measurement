"""Extract frames and space-time data from the high-speed camera videos (.cine) of the plug and slug runs.

Adapted from old_notebooks/save_camera_3h_final.ipynb (crop strips per recording day, gamma display).
Differences: the frame rate and time stamps come from each .cine header (the old notebook assumed
1000 fps; the videos are at 5000 or 8000 fps), frames are kept at full resolution as lossless PNG,
the contrast is normalized once per run (no flicker between frames), and space-time arrays are
computed from every frame.

Output (unpublished experimental data, ignored by Git): data/bronze_layer/camera/run_<id>_<date_identifier>/
    metadata.json           run, video and crop information
    frames/                 cropped frames every --step-ms, f<index>_t<time ms>.png (8-bit, display contrast)
    spacetime.npz           from every frame: t_s, col_mean (frames x columns), probes (frames x rows, one per
                            probe column), probe_cols, mean_img, std_img (raw intensity, cropped)
    contact_sheet.png       frames evenly spaced over the video, with times
    probe_mid.png           height (rows) vs time at the middle probe column
    spacetime.png           column-mean intensity, time (rows) vs position (columns)

Usage:
    python scripts/extract_camera_frames.py --runs 0305_06 0323_02      # date identifiers or run ids
    python scripts/extract_camera_frames.py --all                        # the 58 plug and slug runs
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pims
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RUNS_TABLE = ROOT / "data/gold_layer/objective_criteria_11/processed_v2/objective_flow_pattern.pkl"
OUT_ROOT = ROOT / "data/bronze_layer/camera"

# Crop strips [top, bottom] in pixels per recording day (old_notebooks/save_camera_3h_final.ipynb)
CROP_STRIPS = {"0227": (140, 150), "0228": (95, 190), "0304": (140, 90), "0305": (110, 150),
               "0306": (150, 120), "0312": (155, 120), "0313": (145, 140), "0314": (170, 125),
               "0316": (155, 135), "0318": (140, 140), "0319": (150, 150), "0321": (150, 150),
               "0322": (160, 110), "0323": (175, 90), "0325": (145, 130), "0326": (185, 105)}
CROP_EXCEPTIONS = {("0322", "v2"): (190, 85), ("0325", "v3"): (165, 120)}
DEFAULT_STRIP = (120, 120)

GAMMA = 2.222                  # display gamma (old notebook)
CONTRAST_PERCENTILES = (0.5, 99.5)
PROBE_FRACTIONS = (0.25, 0.50, 0.75)   # probe columns as fractions of the frame width


def cine_name(date_identifier):
    """0305_06 -> 0305_06.cine; 0322_v2_04 -> 0322_04_v2.cine (file naming on the drive)."""
    parts = date_identifier.split("_")
    return f"{date_identifier}.cine" if len(parts) == 2 else f"{parts[0]}_{parts[2]}_{parts[1]}.cine"


def crop_strip(date_identifier):
    parts = date_identifier.split("_")
    day, version = parts[0], (parts[1] if len(parts) == 3 else "")
    return CROP_EXCEPTIONS.get((day, version), CROP_STRIPS.get(day, DEFAULT_STRIP))


def to_display(raw, lo, hi):
    x = np.clip((raw.astype(np.float32) - lo) / (hi - lo), 0, 1)
    return (np.power(x, 1 / GAMMA) * 255).astype(np.uint8)


class RetryVideo:
    """pims Cine reader that reopens the file and retries after transient read errors of the external drive."""

    def __init__(self, path, retries=3):
        self.path, self.retries = path, retries
        self.video = pims.open(path)

    def __len__(self):
        return len(self.video)

    def __getattr__(self, name):
        return getattr(self.video, name)

    def __getitem__(self, i):
        for attempt in range(self.retries + 1):
            try:
                return np.asarray(self.video[i])
            except OSError:
                if attempt == self.retries:
                    raise
                self.video.close()
                self.video = pims.open(self.path)

    def close(self):
        self.video.close()


def extract(run, drive, step_ms, n_sheet, out_root, overwrite):
    di = run.date_identifier
    out = out_root / f"run_{int(run.run_id):03d}_{di}"
    if out.exists() and not overwrite:
        print(f"{di}: exists, skipped ({out})")
        return
    (out / "frames").mkdir(parents=True, exist_ok=True)

    video = RetryVideo(str(Path(drive) / cine_name(di)))
    n, fps = len(video), float(video.frame_rate)
    top, bottom = crop_strip(di)
    first = np.asarray(video[0])            # (rows, columns); video.frame_shape reports them swapped
    rows = slice(top, first.shape[0] - bottom)
    h, w = first[rows, :].shape
    probe_cols = [int(round(f * (w - 1))) for f in PROBE_FRACTIONS]

    # time stamps of each frame relative to the first one (from the header; fall back to index / fps)
    try:                                    # stamps are (datetime to the second, fraction of the second)
        stamps = video.frame_time_stamps
        t = np.array([(d - stamps[0][0]).total_seconds() + (frac - stamps[0][1]) for d, frac in stamps])
        if len(t) != n or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
            raise ValueError
        time_source = "frame_time_stamps"
    except Exception:
        t = np.arange(n) / fps
        time_source = "index/frame_rate"

    # per-run contrast from a sample of frames
    sample = np.stack([np.asarray(video[i])[rows, :] for i in np.linspace(0, n - 1, min(n, 40)).astype(int)])
    lo, hi = (float(v) for v in np.percentile(sample, CONTRAST_PERCENTILES))

    step = max(1, int(round(step_ms * 1e-3 * fps)))
    col_mean = np.empty((n, w), np.float32)
    probes = np.empty((len(probe_cols), n, h), np.uint16)
    s1 = np.zeros((h, w), np.float64)
    s2 = np.zeros((h, w), np.float64)
    saved = []
    for i in range(n):
        crop = np.asarray(video[i])[rows, :]
        col_mean[i] = crop.mean(axis=0)
        for k, c in enumerate(probe_cols):
            probes[k, i] = crop[:, c]
        f = crop.astype(np.float64)
        s1 += f
        s2 += f * f
        if i % step == 0:
            name = f"f{i:05d}_t{t[i] * 1e3:08.3f}ms.png"
            Image.fromarray(to_display(crop, lo, hi)).save(out / "frames" / name, optimize=True)
            saved.append(name)
    video.close()
    mean_img = s1 / n
    std_img = np.sqrt(np.maximum(s2 / n - mean_img ** 2, 0))

    np.savez_compressed(out / "spacetime.npz", t_s=t.astype(np.float64), col_mean=col_mean.astype(np.float16),
                        probes=probes, probe_cols=np.array(probe_cols), mean_img=mean_img.astype(np.float32),
                        std_img=std_img.astype(np.float32))

    # overview images
    idx = np.linspace(0, n - 1, n_sheet).astype(int)
    video = RetryVideo(str(Path(drive) / cine_name(di)))
    fig, axs = plt.subplots(n_sheet, 1, figsize=(8, n_sheet * (8 * h / w + 0.25) + 0.3), squeeze=False,
                            layout="constrained")
    for ax, i in zip(axs[:, 0], idx):
        ax.imshow(to_display(np.asarray(video[i])[rows, :], lo, hi), cmap="gray", vmin=0, vmax=255, aspect="equal")
        ax.set_axis_off()
        ax.set_title(f"t = {t[i] * 1e3:.1f} ms", fontsize=7, pad=1, loc="left")
    video.close()
    fig.suptitle(f"run {int(run.run_id)} ({di}), {run.reduced_pattern}, $j_l$ = {run.jl:.2f} m/s, "
                 f"$j_g$ = {run.jg:.2f} m/s, {fps:.0f} fps", fontsize=8)
    fig.savefig(out / "contact_sheet.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    mid = probes[len(probe_cols) // 2].T  # rows x frames
    fig, ax = plt.subplots(figsize=(8, 2.2))
    ax.imshow(to_display(mid, lo, hi), cmap="gray", aspect="auto", extent=[t[0] * 1e3, t[-1] * 1e3, h, 0])
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("row [px]")
    ax.set_title(f"run {int(run.run_id)} ({di}): column {probe_cols[len(probe_cols) // 2]} vs time", fontsize=8)
    fig.savefig(out / "probe_mid.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 6))
    ax.imshow(col_mean, cmap="gray", aspect="auto", extent=[0, w, t[-1] * 1e3, t[0] * 1e3])
    ax.set_xlabel("column [px]")
    ax.set_ylabel("t [ms]")
    ax.set_title(f"run {int(run.run_id)} ({di}): column mean", fontsize=8)
    fig.savefig(out / "spacetime.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    meta = dict(run_id=int(run.run_id), date_identifier=di, cine=cine_name(di), pattern=run.reduced_pattern,
                jl_m_s=float(run.jl), jg_m_s=float(run.jg), n_frames=n, fps=fps, duration_s=float(t[-1] - t[0]),
                time_source=time_source, frame_shape_raw=list(first.shape), crop_rows=[top, first.shape[0] - bottom],
                crop_shape=[h, w], step_frames=step, step_ms=step / fps * 1e3, n_frames_saved=len(saved),
                display=dict(lo=lo, hi=hi, gamma=GAMMA, percentiles=list(CONTRAST_PERCENTILES)),
                probe_cols=probe_cols, source_script="scripts/extract_camera_frames.py")
    (out / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"{di}: {n} frames at {fps:.0f} fps ({meta['duration_s'] * 1e3:.0f} ms), {len(saved)} saved, crop {h}x{w} -> {out}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--runs", nargs="+", help="date identifiers (e.g. 0305_06) or run ids")
    g.add_argument("--all", action="store_true", help="all plug and slug runs")
    p.add_argument("--drive", default="D:/", help="folder with the .cine files")
    p.add_argument("--step-ms", type=float, default=2.0, help="time between saved frames [ms]")
    p.add_argument("--sheet", type=int, default=8, help="frames in the contact sheet")
    p.add_argument("--out", default=str(OUT_ROOT))
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args(argv)

    warnings.filterwarnings("ignore", category=UserWarning)
    runs = pd.read_pickle(RUNS_TABLE)
    runs = runs[runs.reduced_pattern.isin(["plug", "slug"])]
    if not a.all:
        keys = set(a.runs)
        runs = runs[runs.date_identifier.isin(keys) | runs.run_id.astype(str).isin(keys)]
        missing = keys - set(runs.date_identifier) - set(runs.run_id.astype(str))
        if missing:
            sys.exit(f"not plug/slug runs or unknown: {sorted(missing)}")
    out_root = Path(a.out)
    out_root.mkdir(parents=True, exist_ok=True)
    for run in runs.itertuples():
        extract(run, a.drive, a.step_ms, a.sheet, out_root, a.overwrite)


if __name__ == "__main__":
    main()
