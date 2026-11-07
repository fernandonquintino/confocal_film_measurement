"""
Slug units from the confocal signal of plug and slug flows.

A slug unit is a liquid slug followed by an elongated bubble. At the bottom-wall probe:

* the elongated bubble is a film interval: the detected interface stays close to the wall
  (film state) or is below the detection limit;
* the liquid slug is a cluster of tall readings (second thickness >= h_tall_rel * D_h): liquid
  bridges in unaerated slugs, or the bottoms of the dispersed bubbles of aerated slugs, which
  leave no liquid bridge (Fluids paper, Stage 7; the camera shows aerated slugs filling the
  channel while the probe reads 2-4 mm). Tall readings separated by less than the shortest
  elongated bubble (a length, min_bubble_len_rel * D_h, converted to time with j_m) belong to
  the same slug.

Samples removed as mixer waves (processed_v2) are counted as tall readings: in slug flow the
short tall readings inside the film are mostly the passage of aerated slugs.

Durations are measured directly. Lengths need a translational velocity; here the mixture
superficial velocity j_m = j_l + j_g is used as an approximation (L = j_m t).
"""

import numpy as np
import pandas as pd

from src.flow_pattern_criteria import CHANNEL_SIDE, film_estimate, merge_params

SLUG_PARAMS = {
    'h_tall_rel': 0.25,          # tall reading: t2 / D_h at or above this (same limit as the film state) [-]
    'min_bubble_len_rel': 2.0,   # shortest elongated bubble, L_B,min / D_h; converted to time with j_m [-]
    'min_gap_samples': 2,        # the bubble time limit is never shorter than this many samples [-]
    'min_slug': 0.005,           # slugs shorter than this are isolated readings, not slugs [s]
}

# sensitivity band: all combinations of these values
SLUG_PARAM_GRID = {
    'min_bubble_len_rel': [1.0, 2.0, 4.0],
    'min_slug': [0.0025, 0.005, 0.010],   # 0 would count every isolated tall reading as a slug
}

UNIT_COLUMNS = ['s_start', 's_end', 'b_end', 't_s', 't_b']


def detect_slug_units(df, jm, params=None):
    """
    Slug units of one processed run (columns time, t2 [mm], optional mixer_wave); jm = j_l + j_g [m/s].

    Tall readings separated by a film interval shorter than the shortest elongated bubble,
    max(L_B,min / j_m, min_gap_samples * dt), belong to the same liquid slug.

    Returns (units, info):
        units  DataFrame, one row per complete unit: slug start/end, bubble end [s],
               slug duration t_s and bubble duration t_b [s]
        info   dict: number of slugs, sample interval, record length, bubble time limit
    """

    p = merge_params(SLUG_PARAMS, params, 'SLUG_PARAMS')

    t = df['time'].to_numpy(dtype=float)
    t2 = df['t2'].to_numpy(dtype=float)
    mixer = df['mixer_wave'].to_numpy(dtype=bool) if 'mixer_wave' in df else np.zeros(len(df), bool)

    dt = float(np.nanmedian(np.diff(t)))
    duration = float(t[-1] - t[0] + dt)
    t_bubble_min = max(p['min_bubble_len_rel'] * CHANNEL_SIDE * 1e-3 / jm, p['min_gap_samples'] * dt)
    info = {'dt': dt, 'duration': duration, 'n_slugs': 0, 't_bubble_min': t_bubble_min}
    empty = pd.DataFrame(columns=UNIT_COLUMNS, dtype=float)

    tall = (np.nan_to_num(t2, nan=-1.0) >= p['h_tall_rel'] * CHANNEL_SIDE) | mixer
    idx = np.flatnonzero(tall)

    if len(idx) == 0:
        return empty, info

    # tall readings separated by less than the shortest bubble form one slug
    breaks = np.flatnonzero(np.diff(t[idx]) - dt >= t_bubble_min)
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    slugs = [[t[s], t[e] + dt] for s, e in zip(starts, ends) if t[e] + dt - t[s] >= p['min_slug']]

    if not slugs:
        return empty, info

    slugs = np.array(slugs)
    info['n_slugs'] = len(slugs)

    if len(slugs) < 2:
        return empty, info

    units = pd.DataFrame({'s_start': slugs[:-1, 0], 's_end': slugs[:-1, 1], 'b_end': slugs[1:, 0]})
    units['t_s'] = units['s_end'] - units['s_start']
    units['t_b'] = units['b_end'] - units['s_end']

    return units, info


def unit_statistics(units, info, jl, jg):
    """
    Run-level slug quantities.

    Measured (durations only):
        f_slug        slugs per second of record [Hz]
        t_s_med, t_b_med   median slug and bubble durations [s]
        intermittency fraction of the unit time occupied by the bubble, sum(t_b) / sum(t_s + t_b)
        ratio_bs      median of t_b / t_s per unit = L_B / L_S if nose and tail travel at the same velocity
    Estimated with U = j_m:
        Ls_Dh, Lb_Dh  median slug and bubble lengths relative to D_h
        St            f_slug D_h / j_m
    """

    jm = jl + jg
    d_h = CHANNEL_SIDE * 1e-3
    out = {'n_units': len(units), 'n_slugs': info['n_slugs'], 'duration_s': info['duration'],
           'f_slug': info['n_slugs'] / info['duration'], 'jm': jm}

    if len(units):
        out.update({
            't_s_med': units['t_s'].median(), 't_b_med': units['t_b'].median(),
            'intermittency': units['t_b'].sum() / (units['t_s'] + units['t_b']).sum(),
            'ratio_bs': (units['t_b'] / units['t_s']).median(),
            'Ls_Dh': jm * units['t_s'].median() / d_h, 'Lb_Dh': jm * units['t_b'].median() / d_h,
        })
    else:
        out.update(dict.fromkeys(['t_s_med', 't_b_med', 'intermittency', 'ratio_bs', 'Ls_Dh', 'Lb_Dh'], np.nan))

    out['St'] = out['f_slug'] * d_h / jm
    return out


def film_below_bubbles(df, units, n_bins=20, min_samples=5):
    """
    Bottom film below the elongated bubbles as a function of the position in the bubble,
    x = (t - s_end) / t_b (0: just behind the slug, 1: next slug).

    Each bubble with at least min_samples film values contributes its bin means; returns the
    bin centres, the median over bubbles of the bin means [mm] and the number of bubbles per bin.
    Uses the film estimate (detected interface, or film merged with the wall reading).
    """

    t = df['time'].to_numpy(dtype=float)
    film = film_estimate(df)
    edges = np.linspace(0, 1, n_bins + 1)
    rows = []

    for u in units.itertuples():
        m = (t >= u.s_end) & (t < u.b_end) & np.isfinite(film)
        if m.sum() < min_samples:
            continue
        x = (t[m] - u.s_end) / u.t_b
        b = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
        means = np.full(n_bins, np.nan)
        for k in np.unique(b):
            means[k] = np.mean(film[m][b == k])
        rows.append(means)

    centres = 0.5 * (edges[1:] + edges[:-1])

    if not rows:
        return centres, np.full(n_bins, np.nan), np.zeros(n_bins, int)

    rows = np.array(rows)
    return centres, np.nanmedian(rows, axis=0), np.sum(np.isfinite(rows), axis=0)
