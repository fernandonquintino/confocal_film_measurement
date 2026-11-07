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
    contact_sheet.pdf       saved frames evenly spaced over the video, with camera times
    probe_mid.pdf           height (rows) vs camera time at the middle probe column
    spacetime.pdf           column-mean intensity, camera time (rows) vs position (columns)

Camera time t_cam starts at the first recorded frame of each video. The camera and the confocal sensor
were not synchronized, so t_cam = 0 and the confocal t = 0 are different instants. Figures identify the
runs by pattern and superficial velocities; the run id is kept in the folder name for traceability.

Usage:
    python scripts/extract_camera_frames.py --runs 0305_06 0323_02      # date identifiers or run ids
    python scripts/extract_camera_frames.py --all                        # the 58 plug and slug runs
    python scripts/extract_camera_frames.py --all --overviews-only       # rebuild the overview PDFs only
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
sys.path.insert(0, str(ROOT))
from src.flow_pattern_report import (FIG_SIZE, FIG_SIZE_WIDE, T_CAMERA_LABEL, TIME_ORIGIN_NOTE,   # noqa: E402
                                     run_label)

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


def overviews(out, n_sheet=8):
    """Overview PDFs of one extracted run, built from the saved frames and spacetime.npz (no video needed).

    Times are camera times: zero at the first recorded frame, not synchronized with the confocal record.
    Old PNG overviews of earlier versions of this script are removed.
    """
    meta = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    z = np.load(out / "spacetime.npz")
    t_ms = z["t_s"] * 1e3
    lo, hi = meta["display"]["lo"], meta["display"]["hi"]
    h, w = meta["crop_shape"]
    title = f"{run_label(meta['pattern'], meta['jl_m_s'], meta['jg_m_s'])}, {meta['fps']:.0f} fps"
    note = dict(t=TIME_ORIGIN_NOTE, fontsize=6, color="0.35")

    frames = sorted((out / "frames").glob("f*_t*ms.png"))
    pick = [frames[i] for i in np.linspace(0, len(frames) - 1, min(n_sheet, len(frames))).astype(int)]
    fig, axs = plt.subplots(len(pick), 1, figsize=(FIG_SIZE_WIDE[0], len(pick) * (FIG_SIZE_WIDE[0] * h / w + 0.22) + 0.5),
                            squeeze=False, layout="constrained")
    for ax, f in zip(axs[:, 0], pick):
        ax.imshow(np.asarray(Image.open(f)), cmap="gray", vmin=0, vmax=255, aspect="equal")
        ax.set_axis_off()
        ax.set_title(f"$t_{{cam}}$ = {float(f.stem.split('_t')[1][:-2]):.1f} ms", fontsize=7, pad=1, loc="left")
    fig.suptitle(title, fontsize=8)
    fig.get_layout_engine().set(h_pad=0.02)
    fig.supxlabel(**note)
    fig.savefig(out / "contact_sheet.pdf", dpi=200)
    plt.close(fig)

    mid = len(z["probe_cols"]) // 2
    fig, ax = plt.subplots(figsize=(FIG_SIZE_WIDE[0], 2.4), layout="constrained")
    ax.imshow(to_display(z["probes"][mid].T, lo, hi), cmap="gray", aspect="auto", extent=[t_ms[0], t_ms[-1], h, 0])
    ax.set_xlabel(T_CAMERA_LABEL)
    ax.set_ylabel("row [px]")
    ax.set_title(f"{title}: column {int(z['probe_cols'][mid])} px vs time", fontsize=8, loc="left")
    fig.supxlabel(**note)
    fig.savefig(out / "probe_mid.pdf", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(FIG_SIZE[0], 1.4 * FIG_SIZE[1]), layout="constrained")
    ax.imshow(z["col_mean"].astype(np.float32), cmap="gray", aspect="auto", extent=[0, w, t_ms[-1], t_ms[0]])
    ax.set_xlabel("column [px]")
    ax.set_ylabel(T_CAMERA_LABEL)
    ax.set_title(f"{title}\ncolumn-mean intensity", fontsize=7, loc="left")
    fig.supxlabel(**note)
    fig.savefig(out / "spacetime.pdf", dpi=200)
    plt.close(fig)

    for old in ("contact_sheet.png", "probe_mid.png", "spacetime.png"):
        (out / old).unlink(missing_ok=True)


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

    meta = dict(run_id=int(run.run_id), date_identifier=di, cine=cine_name(di), pattern=run.reduced_pattern,
                jl_m_s=float(run.jl), jg_m_s=float(run.jg), n_frames=n, fps=fps, duration_s=float(t[-1] - t[0]),
                time_source=time_source, frame_shape_raw=list(first.shape), crop_rows=[top, first.shape[0] - bottom],
                crop_shape=[h, w], step_frames=step, step_ms=step / fps * 1e3, n_frames_saved=len(saved),
                display=dict(lo=lo, hi=hi, gamma=GAMMA, percentiles=list(CONTRAST_PERCENTILES)),
                probe_cols=probe_cols, time_origin="first recorded frame (camera clock); not synchronized with the "
                "confocal record", source_script="scripts/extract_camera_frames.py")
    (out / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    overviews(out, n_sheet)
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
    p.add_argument("--overviews-only", action="store_true",
                   help="only rebuild the overview PDFs of already extracted runs (no video needed)")
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
        if a.overviews_only:
            folder = out_root / f"run_{int(run.run_id):03d}_{run.date_identifier}"
            if (folder / "spacetime.npz").exists():
                overviews(folder, a.sheet)
                print(f"{run.date_identifier}: overviews -> {folder}")
            continue
        extract(run, a.drive, a.step_ms, a.sheet, out_root, a.overwrite)


if __name__ == "__main__":
    main()
