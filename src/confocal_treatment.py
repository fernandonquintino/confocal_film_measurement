import numpy as np
import pandas as pd

from src.flow_pattern_criteria import CHANNEL_SIDE, merge_params

# =========================================================
# Mixer wave removal (processed_v2)
# =========================================================
#
# In separated flows (stratified, wavy, annular) the liquid accumulated at
# the mixer is released as short waves that almost fill the channel. At the
# probe they appear as isolated tall readings (t2 of 6 mm to 8 mm) lasting
# 1 ms to 3 ms inside a thin film. Real liquid slugs/plugs last 15 ms to
# 130 ms, and the liquid of bubbly flow is surrounded by other tall readings.
#
# An event is removed only when it is
#   * short    : longest contiguous tall stretch < MAX_DURATION
#   * in a film: median of the valid t2 around it < CONTEXT_MAX_FILM
#   * isolated : fraction of tall readings around it < CONTEXT_MAX_TALL
#
# At 100 Hz one sample lasts 10 ms, so the duration test cannot be met and
# those runs are left untouched (spurs and short slugs are not resolvable).

# Lengths are dimensionless (``_rel`` = length / CHANNEL_SIDE, the nominal inner
# side H = D_h of the square channel; value in mm in the comment). Times stay
# dimensional: they are tied to the sample interval.
MIXER_WAVE_PARAMS = {
    'h_tall_rel': 0.5,              # tall reading: t2 / H >= h_tall_rel (4.0 mm) [-]
    'h_film_rel': 0.25,             # neighbouring samples with t2 / H >= h_film_rel belong to the event (2.0 mm) [-]
    'max_gap': 0.02,                # tall readings closer than this form one event [s]
    'max_duration': 0.005,          # longest contiguous tall stretch of a mixer wave [s]
    'context_window': 0.2,          # window before/after the event [s]
    'context_max_film_rel': 0.25,   # median valid t2 / H around the event (2.0 mm) [-]
    'context_max_tall': 0.1,        # fraction of tall readings around the event [-]
}


def detect_mixer_waves(df, params=None):
    """
    Boolean mask (aligned with ``df``) of the samples that belong to mixer
    waves. See the module header for the criteria.
    """

    p = merge_params(MIXER_WAVE_PARAMS, params, 'MIXER_WAVE_PARAMS')

    # Length parameters in mm (exact: CHANNEL_SIDE is a power of two)
    p['h_tall'] = p['h_tall_rel'] * CHANNEL_SIDE
    p['h_film'] = p['h_film_rel'] * CHANNEL_SIDE
    p['context_max_film'] = p['context_max_film_rel'] * CHANNEL_SIDE

    time = df['time'].to_numpy(dtype=float)
    t2 = df['t2'].to_numpy(dtype=float)

    mask = np.zeros(len(df), dtype=bool)

    if len(df) < 2:
        return pd.Series(mask, index=df.index)

    dt = np.nanmedian(np.diff(time))

    tall = np.nan_to_num(t2, nan=-1) >= p['h_tall']
    tall_idx = np.where(tall)[0]

    if len(tall_idx) == 0:
        return pd.Series(mask, index=df.index)

    # Group tall readings into events
    breaks = np.where(np.diff(time[tall_idx]) > p['max_gap'])[0] + 1

    for group in np.split(tall_idx, breaks):

        start, end = group[0], group[-1]

        # Longest contiguous tall stretch inside the event
        edges = np.diff(np.concatenate([[0], tall[start:end + 1].astype(int), [0]]))
        longest = (np.where(edges == -1)[0] - np.where(edges == 1)[0]).max() * dt

        if longest >= p['max_duration']:
            continue

        # Context: valid t2 before and after the event
        context = (
            ((time >= time[start] - p['context_window']) & (time < time[start]))
            | ((time > time[end]) & (time <= time[end] + p['context_window']))
        )
        context_t2 = t2[context]
        context_t2 = context_t2[~np.isnan(context_t2)]

        if len(context_t2) == 0:
            continue

        if np.median(context_t2) >= p['context_max_film']:
            continue

        if np.mean(context_t2 >= p['h_tall']) >= p['context_max_tall']:
            continue

        # Extend the event to the neighbouring non-film samples (rise and fall)
        lo, hi = start, end

        while lo > 0 and np.nan_to_num(t2[lo - 1], nan=-1) >= p['h_film']:
            lo -= 1

        while hi < len(t2) - 1 and np.nan_to_num(t2[hi + 1], nan=-1) >= p['h_film']:
            hi += 1

        mask[lo:hi + 1] = True

    return pd.Series(mask, index=df.index)


def remove_mixer_waves(df, params=None):
    """
    Copy of ``df`` with the mixer-wave samples removed from t2 (set to NaN)
    and flagged in a boolean ``mixer_wave`` column, so the features can
    exclude them instead of counting them as undetected.
    """

    df = df.copy()

    mask = detect_mixer_waves(df, params)

    df.loc[mask, 't2'] = np.nan
    df['mixer_wave'] = mask

    return df
