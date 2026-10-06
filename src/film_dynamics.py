import math

import numpy as np
import pandas as pd

from scipy.signal import find_peaks, lombscargle, peak_widths

from src.flow_pattern_criteria import FEATURE_PARAMS, film_estimate, load_run

# =========================================================
# Film dynamics (churn vs. wavy)
# =========================================================
#
# All features are computed on the bottom film estimate (detected t2 below
# h_liq, or the film merged with the wall reading), averaged into 10 ms
# bins. The common 100 Hz grid makes the 100 Hz and kHz runs comparable.
# Gaps up to MAX_FILL_BINS are interpolated; longer gaps split the signal
# into segments.

BIN_S = 0.01
MAX_FILL_BINS = 3

SPEC_FREQS = np.linspace(0.5, 50.0, 400)   # Hz (50 Hz = Nyquist of the 10 ms grid)
SPEC_LOW_HZ = 5.0
SPEC_MIN_BINS = 200

WAVE_BASE_WINDOW = 50       # bins (0.5 s) for the rolling base film
WAVE_BASE_QUANTILE = 0.10
WAVE_MIN_PROMINENCE = 0.1   # mm
WAVE_REL_PROMINENCE = 0.3   # x base film
WAVE_MIN_BASE = 0.05        # mm, floor of the base film in relative heights
WAVE_MIN_SEGMENT = 50       # bins

PERM_ORDER = 4
PERM_MIN_PATTERNS = 100

WINDOW_MIN_VALID = 0.5
WINDOW_MIN_COUNT = 4

DYNAMICS_FEATURES = [
    'spec_dominant_hz', 'spec_centroid_hz', 'spec_entropy', 'spec_flatness', 'spec_low_share',
    'wave_freq_hz', 'wave_height_rel', 'wave_large_frac', 'wave_interval_cv', 'wave_asymmetry',
    'perm_entropy', 'window_cv_05s', 'window_cv_1s', 'gap_burstiness', 'gap_fano_05s',
]

# =========================================================
# Signal preparation
# =========================================================

def _drop_mixer_waves(df):

    if 'mixer_wave' in df.columns:
        df = df.loc[~df['mixer_wave'].astype(bool)].reset_index(drop=True)

    return df


def resample_film(df, h_liq=FEATURE_PARAMS['h_liq'], bin_s=BIN_S):
    """
    Bottom film estimate averaged into ``bin_s`` bins. Returns the bin times
    and values (NaN where no valid sample fell in the bin).
    """

    df = _drop_mixer_waves(df)

    estimate = film_estimate(df, h_liq)
    bins = np.floor(df['time'].to_numpy() / bin_s).astype(int)

    valid = ~np.isnan(estimate)
    n_bins = bins.max() + 1 if len(bins) else 0

    sums = np.bincount(bins[valid], weights=estimate[valid], minlength=n_bins)
    counts = np.bincount(bins[valid], minlength=n_bins)

    values = np.full(n_bins, np.nan)
    values[counts > 0] = sums[counts > 0] / counts[counts > 0]

    return np.arange(n_bins) * bin_s, values


def _fill_short_gaps(values, max_bins=MAX_FILL_BINS):

    s = pd.Series(values)
    filled = s.interpolate(limit=max_bins, limit_area='inside')

    # interpolate(limit=...) also fills the first bins of longer gaps: undo it
    is_nan = s.isna().to_numpy()
    edges = np.diff(np.concatenate([[0], is_nan.astype(int), [0]]))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0]

    out = filled.to_numpy(copy=True)

    for start, end in zip(starts, ends):
        if end - start > max_bins:
            out[start:end] = np.nan

    return out


def _segments(values, min_length):

    valid = ~np.isnan(values)
    edges = np.diff(np.concatenate([[0], valid.astype(int), [0]]))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0]

    return [(s, e) for s, e in zip(starts, ends) if e - s >= min_length]

# =========================================================
# a) Spectral features
# =========================================================

def spectral_features(times, values):
    """
    Lomb-Scargle periodogram of the film (valid bins only, 0.5 Hz to 50 Hz),
    normalised to unit power.
    """

    valid = ~np.isnan(values)

    keys = ['spec_dominant_hz', 'spec_centroid_hz', 'spec_entropy', 'spec_flatness', 'spec_low_share']

    if valid.sum() < SPEC_MIN_BINS or np.nanstd(values) == 0:
        return dict.fromkeys(keys, np.nan)

    t = times[valid]
    y = values[valid] - values[valid].mean()

    power = lombscargle(t, y, 2 * np.pi * SPEC_FREQS)
    power = np.clip(power, 1e-30, None)
    p = power / power.sum()

    return {
        'spec_dominant_hz': SPEC_FREQS[np.argmax(power)],
        'spec_centroid_hz': float(np.sum(SPEC_FREQS * p)),
        'spec_entropy': float(-np.sum(p * np.log(p)) / np.log(len(p))),
        'spec_flatness': float(np.exp(np.mean(np.log(power))) / np.mean(power)),
        'spec_low_share': float(p[SPEC_FREQS < SPEC_LOW_HZ].sum()),
    }

# =========================================================
# b) Wave events
# =========================================================

def wave_features(values, bin_s=BIN_S):
    """
    Individual waves above a rolling base film (10th percentile over 0.5 s)
    in the gap-filled segments.
    """

    keys = ['wave_freq_hz', 'wave_height_rel', 'wave_large_frac', 'wave_interval_cv', 'wave_asymmetry']

    filled = _fill_short_gaps(values)
    segments = _segments(filled, WAVE_MIN_SEGMENT)

    if not segments:
        return dict.fromkeys(keys, np.nan)

    heights_rel, intervals, rises, falls = [], [], [], []
    n_waves, total_time = 0, 0.0

    for start, end in segments:

        y = filled[start:end]

        base = (
            pd.Series(y)
            .rolling(WAVE_BASE_WINDOW, center=True, min_periods=10)
            .quantile(WAVE_BASE_QUANTILE)
            .to_numpy()
        )
        base_floor = np.maximum(base, WAVE_MIN_BASE)

        prominence = max(WAVE_MIN_PROMINENCE, WAVE_REL_PROMINENCE * np.nanmedian(base_floor))
        peaks, _ = find_peaks(y, prominence=prominence)

        total_time += (end - start) * bin_s

        if len(peaks) == 0:
            continue

        n_waves += len(peaks)

        heights_rel.extend((y[peaks] - base[peaks]) / base_floor[peaks])
        intervals.extend(np.diff(peaks) * bin_s)

        _, _, left, right = peak_widths(y, peaks, rel_height=0.5)
        rises.extend(peaks - left)
        falls.extend(right - peaks)

    heights_rel = np.asarray(heights_rel, dtype=float)
    intervals = np.asarray(intervals, dtype=float)

    return {
        'wave_freq_hz': n_waves / total_time if total_time > 0 else np.nan,
        'wave_height_rel': float(np.nanmedian(heights_rel)) if len(heights_rel) else np.nan,
        'wave_large_frac': float(np.nanmean(heights_rel > 2)) if len(heights_rel) else np.nan,
        'wave_interval_cv': float(intervals.std() / intervals.mean()) if len(intervals) >= 5 else np.nan,
        'wave_asymmetry': float(np.median(rises) / np.median(falls)) if len(rises) >= 5 and np.median(falls) > 0 else np.nan,
    }

# =========================================================
# c) Complexity and unsteadiness
# =========================================================

def permutation_entropy(values, order=PERM_ORDER):
    """
    Normalised permutation entropy, pooled over the gap-filled segments
    (0 = fully regular, 1 = fully random ordering).
    """

    filled = _fill_short_gaps(values)
    counts = {}

    for start, end in _segments(filled, order + 1):

        windows = np.lib.stride_tricks.sliding_window_view(filled[start:end], order)
        codes = np.argsort(windows, axis=1, kind='stable') @ (order ** np.arange(order))

        for code, n in zip(*np.unique(codes, return_counts=True)):
            counts[code] = counts.get(code, 0) + n

    total = sum(counts.values())

    if total < PERM_MIN_PATTERNS:
        return np.nan

    p = np.array(list(counts.values())) / total

    return float(-np.sum(p * np.log(p)) / np.log(math.factorial(order)))


def window_cv(values, window_bins):
    """
    Coefficient of variation of the window means (windows with at least
    WINDOW_MIN_VALID valid bins): slow alternation of thick and thin film.
    """

    n_windows = len(values) // window_bins

    means = []

    for i in range(n_windows):

        w = values[i * window_bins:(i + 1) * window_bins]

        if np.mean(~np.isnan(w)) >= WINDOW_MIN_VALID:
            means.append(np.nanmean(w))

    if len(means) < WINDOW_MIN_COUNT:
        return np.nan

    means = np.asarray(means)

    return float(means.std() / means.mean()) if means.mean() > 0 else np.nan


def gap_burstiness(df, window_s=0.5):
    """
    Burstiness B = (sigma - mu) / (sigma + mu) of the intervals between the
    starts of undetected t2 gaps (-1 periodic, 0 random, -> 1 clustered) and
    the Fano factor of the gap counts in ``window_s`` windows.
    """

    df = _drop_mixer_waves(df)

    time = df['time'].to_numpy()
    missing = df['t2'].isna().to_numpy()

    edges = np.diff(np.concatenate([[0], missing.astype(int)]))
    gap_starts = time[np.where(edges == 1)[0]]

    if len(gap_starts) < 10:
        return {'gap_burstiness': np.nan, 'gap_fano_05s': np.nan}

    intervals = np.diff(gap_starts)
    mu, sigma = intervals.mean(), intervals.std()

    counts, _ = np.histogram(gap_starts, bins=np.arange(0, time[-1] + window_s, window_s))

    return {
        'gap_burstiness': float((sigma - mu) / (sigma + mu)) if sigma + mu > 0 else np.nan,
        'gap_fano_05s': float(counts.var() / counts.mean()) if len(counts) >= WINDOW_MIN_COUNT and counts.mean() > 0 else np.nan,
    }

# =========================================================
# All features
# =========================================================

def dynamics_features(df):
    """
    All film dynamics features of one run.
    """

    times, values = resample_film(df)

    return {
        **spectral_features(times, values),
        **wave_features(values),
        'perm_entropy': permutation_entropy(values),
        'window_cv_05s': window_cv(values, int(round(0.5 / BIN_S))),
        'window_cv_1s': window_cv(values, int(round(1.0 / BIN_S))),
        **gap_burstiness(df),
    }


def dynamics_table(run_ids, base_path, trim=False, raw_path=None):
    """
    Film dynamics feature table (one row per run).
    """

    return pd.DataFrame([
        {'run_id': run_id, **dynamics_features(load_run(run_id, base_path, trim=trim, raw_path=raw_path))}
        for run_id in run_ids
    ])
