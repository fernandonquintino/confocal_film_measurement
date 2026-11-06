import numpy as np
import pandas as pd

from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, kurtosis, skew

# =========================================================
# Constants
# =========================================================

# Measured geometry and optics, used by the sensor treatment (dimensional) -----
CHANNEL_HEIGHT = 8.02   # processed t2 cap: inner side measured at the probe [mm]
T1_THRESHOLD = 0.985    # wall reference of the tube-thickness correction [mm]
N_WATER = 1.333
N_BK7 = 1.519

# Nominal geometry, used ONLY to make the empirical criterion lengths dimensionless.
# Square channel: hydraulic diameter D_h = 4A/P = H. Parameters and features
# ending in ``_rel`` are lengths divided by CHANNEL_SIDE (dimensionless); those
# ending in ``_mm`` are in mm. H = 8 (a power of two) makes the conversion
# value_rel * CHANNEL_SIDE exact in floating point, so the classification is
# identical to the former thresholds in mm.
CHANNEL_SIDE = 8.0      # nominal inner side H = D_h [mm]

MIN_KDE_SAMPLES = 50    # minimum detected samples to estimate a KDE
LOW_SAMPLE_LIMIT = 1000 # runs shorter than this are flagged

KDE_GRID_MAX_REL = 1.0625   # KDE grid upper limit / H (8.5 mm) [-]
KDE_MIN_STD_REL = 1.25e-4   # detected t2 with a smaller std / H has no KDE (1e-3 mm) [-]
KDE_GRID = np.linspace(0, KDE_GRID_MAX_REL * CHANNEL_SIDE, 851)   # [mm]

GAP_LONG = 0.1          # undetected gaps longer than this are "long blackouts" [s]

STATES = ['undetected', 'film', 'bubble', 'liquid']

# Feature groups (used by the missing-data and clustering analyses)
MISSING_FEATURES = [
    'f_undetected', 'gap_freq_hz', 'gap_mean_ms', 'gap_long_frac', 'gap_near_interface',
]

DETECTED_FEATURES = [
    'f_film_det', 'f_bubble_det', 'f_liquid_det', 'bubble_share', 'f_wall_merge',
    'bridge_freq_hz', 'bridge_mean_duration_s', 'n_peaks', 'peak_film_mm',
    'film_median_mm', 'film_cv', 'bimodality_coefficient',
    'film_est_mean_mm', 'film_est_std_mm', 'film_est_cv',
    'film_jump_frac', 'film_lag_corr',
]

# Parameters used to compute the run features. Lengths are dimensionless
# (``_rel`` = length / CHANNEL_SIDE; value in mm in the comment); times stay
# dimensional, as they are tied to the sample interval and record length.
FEATURE_PARAMS = {
    'h_liq_rel': 0.9375,           # t2 / H >= h_liq_rel -> liquid bridges the channel (7.5 mm) [-]
    'h_film_rel': 0.25,            # t2 / H <  h_film_rel -> wall film (2.0 mm) [-]
    'min_dwell': 0.005,            # minimum duration of a liquid bridge [s]
    'kde_prominence': 0.05,        # KDE peak prominence relative to its maximum [-]
    'jump_lag': 0.01,              # time lag of the film jump statistics [s]
    'jump_min_change_rel': 0.0375, # film change / H counted as a jump within jump_lag (0.3 mm) [-]
}

# Thresholds of the objective criterion (lengths dimensionless, see FEATURE_PARAMS)
CRITERION_PARAMS = {
    'f_bridge_hz': 0.1,             # liquid bridge frequency -> intermittent [Hz]
    'bubble_share_min': 0.5,        # bubble / (bubble + film) -> bubbly [-]
    'plug_min_duration': 0.02,      # mean liquid bridge duration -> plug [s]
    'plug_min_f_liquid': 0.1,       # time fraction of liquid bridges -> plug [-]
    'annular_max_mean_rel': 0.0125, # thin film: mean bottom film / H -> annular (0.1 mm) [-]
    'annular_max_std_rel': 0.025,   # uniform film: bottom film std / H -> annular (0.2 mm) [-]
    'dispersed_min_peaks': 2,       # KDE peaks of a continuous film -> dispersed [-]
    'dispersed_min_mean_rel': 0.05, # thick film: mean bottom film / H -> dispersed (with jumps) (0.4 mm) [-]
    'dispersed_min_jump': 0.4,      # share of abrupt film changes -> dispersed [-]
    'dispersed_min_detected': 0.15, # film detected in a meaningful share of samples -> dispersed [-]
    'agitated_min_homog_ratio': 1.0, # measured / homogeneous dp/dz: phases strongly mixed [-]
    'dispersed_max_phi2': 5.5,      # liquid-reference multiplier close to the liquid-alone reference -> dispersed [-]
    'liquid_rich_max_phi2': 20.0,   # thin film seen, but multiplier between churn/dispersed (<= 15.1) and thin-film wavy (>= 24.8) -> churn/dispersed [-]
    # Refit in objective_criteria_10 with agitated_min_homog_ratio fixed at 1 (was 5.0 and 15.0)
    'liquid_rich_max_film_rel': 0.0625, # thin film seen by the probe: mean bottom film / H (liquid-rich branch) (0.5 mm) [-]
}

# Length keys in mm of objective_criteria_1..10, replaced by the ``_rel`` keys above
_RENAMED_MM_KEYS = {
    'h_liq': 'h_liq_rel', 'h_film': 'h_film_rel', 'jump_min_change': 'jump_min_change_rel',
    'annular_max_mean': 'annular_max_mean_rel', 'annular_max_std': 'annular_max_std_rel',
    'dispersed_min_mean': 'dispersed_min_mean_rel', 'liquid_rich_max_film': 'liquid_rich_max_film_rel',
}


def merge_params(defaults, params, name):
    """
    ``defaults`` overridden by ``params``. Unknown keys raise a KeyError, so that
    an override with a former length key in mm is never silently ignored.
    """

    params = dict(params or {})
    unknown = set(params) - set(defaults)

    if unknown:
        hints = {k: _RENAMED_MM_KEYS[k] for k in unknown if k in _RENAMED_MM_KEYS}
        raise KeyError(
            f'unknown {name} key(s) {sorted(unknown)}'
            + (f'; lengths are now dimensionless (mm / CHANNEL_SIDE): use {hints}' if hints else '')
        )

    return {**defaults, **params}

# Visual patterns merged into a target class (the probe sees waves in all
# visual stratified runs, so stratified is evaluated as wavy)
REFERENCE_MAP = {'stratified': 'wavy'}

# Runs below this liquid superficial velocity are excluded from the analysis [m/s]
JL_MIN = 0.035

# Kept for backward compatibility (H_LIQ and H_FILM in mm)
H_LIQ = FEATURE_PARAMS['h_liq_rel'] * CHANNEL_SIDE
H_FILM = FEATURE_PARAMS['h_film_rel'] * CHANNEL_SIDE
MIN_DWELL = FEATURE_PARAMS['min_dwell']
KDE_PROMINENCE = FEATURE_PARAMS['kde_prominence']

F_BRIDGE_HZ = CRITERION_PARAMS['f_bridge_hz']
BUBBLE_SHARE_MIN = CRITERION_PARAMS['bubble_share_min']
PLUG_MIN_DURATION = CRITERION_PARAMS['plug_min_duration']
PLUG_MIN_F_LIQUID = CRITERION_PARAMS['plug_min_f_liquid']
DISPERSED_MIN_PEAKS = CRITERION_PARAMS['dispersed_min_peaks']

OBJECTIVE_PATTERNS = ['bubbly', 'plug', 'slug', 'churn', 'dispersed', 'wavy', 'annular']


def reference_pattern(visual):
    """
    Map visual patterns onto the target classes (``REFERENCE_MAP``).
    """

    return pd.Series(visual).replace(REFERENCE_MAP).to_numpy()

# =========================================================
# Loading
# =========================================================

def _trim(df):

    first_valid = df['t1'].first_valid_index()

    if first_valid is not None:
        df = df.loc[first_valid:].reset_index(drop=True)
        df['time'] = df['time'] - df['time'].iloc[0]

    return df


def load_run(run_id, base_path, trim=False, raw_path=None):
    """
    Load a confocal run (processed or raw).

    ``trim=True`` drops the leading samples without a valid t1 and resets
    time to 0 at the first valid t1 (already done in the processed runs).

    ``raw_path`` (for processed runs, whose t1 is clipped at the tube
    threshold) adds the unclipped first thickness of the raw run, trimmed
    the same way, as column ``t1_raw``. It carries the film thinner than the
    detection limit, which merges with the wall reading.

    Expected columns:
        time, t1, t2, t3
    """

    df = pd.read_csv(base_path / f'confocal_run_{run_id}.csv')

    if trim:
        df = _trim(df)

    if raw_path is not None:
        raw = _trim(pd.read_csv(raw_path / f'confocal_run_{run_id}.csv'))

        if len(raw) != len(df):
            raise ValueError(f'run {run_id}: raw ({len(raw)}) and processed ({len(df)}) lengths differ')

        df['t1_raw'] = raw['t1'].to_numpy()

    return df


def load_processed_run(run_id, base_path):
    """
    Load a processed confocal run (time starts at the first valid t1).
    """

    return load_run(run_id, base_path)

# =========================================================
# Sample states
# =========================================================

def sample_states(df, h_liq=H_LIQ, h_film=H_FILM):
    """
    Label each sample with a physical state, based on t2:

    * undetected : t2 missing -> film thinner than the detection limit or
                   interface too steep to reflect (or no signal at all when
                   t1 is also missing, e.g. untrimmed raw runs)
    * film       : t2 <  h_film -> wall film (separated flow or under a bubble)
    * bubble     : h_film <= t2 < h_liq -> bubble or bubble nose/tail
    * liquid     : t2 >= h_liq  -> liquid bridges the channel (liquid slug or
                   continuous liquid of bubbly flow)
    """

    t2 = df['t2']

    states = np.select(
        [t2.isna(), t2 < h_film, t2 < h_liq],
        ['undetected', 'film', 'bubble'],
        default='liquid',
    )

    return pd.Series(states, index=df.index)


def _segments(mask, time):
    """
    Contiguous True segments of a boolean mask -> list of (start, end) indices
    and durations (end time - start time + one sample).
    """

    mask = np.asarray(mask, dtype=bool)
    time = np.asarray(time, dtype=float)

    if not mask.any():
        return [], np.array([])

    edges = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0] - 1

    dt = np.nanmedian(np.diff(time)) if len(time) > 1 else 0.0
    durations = time[ends] - time[starts] + dt

    return list(zip(starts, ends)), durations

# =========================================================
# Features
# =========================================================

def bimodality_coefficient(x):
    """
    Sarle's bimodality coefficient. Values above 5/9 (~0.555) suggest a
    bimodal (or strongly skewed) distribution.
    """

    x = np.asarray(x, dtype=float)
    n = len(x)

    if n < 4 or np.std(x) == 0:
        return np.nan

    g = skew(x)
    k = kurtosis(x)  # excess kurtosis

    return (g ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def kde_peaks(x, grid=KDE_GRID, prominence=KDE_PROMINENCE):
    """
    KDE of the detected film thickness and its significant peaks.

    Returns (density, peak_positions, peak_heights).
    """

    x = np.asarray(x, dtype=float)

    if len(x) < MIN_KDE_SAMPLES or np.std(x) < KDE_MIN_STD_REL * CHANNEL_SIDE:
        return None, np.array([]), np.array([])

    density = gaussian_kde(x, bw_method='scott')(grid)

    peaks, _ = find_peaks(
        np.concatenate([[0], density, [0]]),
        prominence=prominence * density.max(),
    )
    peaks = peaks - 1

    return density, grid[peaks], density[peaks]


def run_features(df, params=None):
    """
    Objective features of one confocal run.

    ``params`` overrides entries of ``FEATURE_PARAMS``. The signals stay in mm;
    the dimensionless length parameters are converted with CHANNEL_SIDE.
    """

    p = merge_params(FEATURE_PARAMS, params, 'FEATURE_PARAMS')

    # Length parameters in mm (exact: CHANNEL_SIDE is a power of two)
    p['h_liq'] = p['h_liq_rel'] * CHANNEL_SIDE
    p['h_film'] = p['h_film_rel'] * CHANNEL_SIDE
    p['jump_min_change'] = p['jump_min_change_rel'] * CHANNEL_SIDE

    # Samples removed as mixer waves (processed_v2) are excluded entirely
    n_mixer_wave = 0

    if 'mixer_wave' in df.columns:
        n_mixer_wave = int(df['mixer_wave'].sum())
        df = df.loc[~df['mixer_wave'].astype(bool)].reset_index(drop=True)

    time = df['time'].to_numpy()
    t2 = df['t2']
    states = sample_states(df, h_liq=p['h_liq'], h_film=p['h_film'])

    n = len(df)
    duration = time[-1] - time[0] if n > 1 else 0.0

    feat = {
        'n_samples': n,
        'duration_s': duration,
        'sampling_hz': 1 / np.nanmedian(np.diff(time)) if n > 1 else np.nan,
        'low_sample_flag': n < LOW_SAMPLE_LIMIT,
        'f_no_signal': float(df['t1'].isna().mean()),
        'n_mixer_wave': n_mixer_wave,
    }

    # -----------------------------------------------------
    # State fractions
    # -----------------------------------------------------

    for state in STATES:
        feat[f'f_{state}'] = float((states == state).mean())

    detected = 1 - feat['f_undetected']

    # Same fractions among detected samples only
    for state in ['film', 'bubble', 'liquid']:
        feat[f'f_{state}_det'] = (
            feat[f'f_{state}'] / detected if detected > 0 else np.nan
        )

    # Share of non-bridged detected interfaces that are bubbles (not wall film)
    bubble_film = feat['f_bubble'] + feat['f_film']
    feat['bubble_share'] = feat['f_bubble'] / bubble_film if bubble_film > 0 else 0.0

    # Thin film merged with the wall reading. Processed t1 is clipped at the
    # threshold, raw t1 exceeds it; both are caught by the >= comparison.
    feat['f_wall_merge'] = float((df['t1'] >= T1_THRESHOLD - 1e-9).mean())

    # -----------------------------------------------------
    # Intermittency (liquid bridges)
    # -----------------------------------------------------

    _, bridge_dur = _segments(states == 'liquid', time)
    bridge_dur = bridge_dur[bridge_dur >= p['min_dwell']]

    feat['bridge_count'] = len(bridge_dur)
    feat['bridge_freq_hz'] = len(bridge_dur) / duration if duration > 0 else np.nan
    feat['bridge_mean_duration_s'] = bridge_dur.mean() if len(bridge_dur) else 0.0
    feat['bridge_duty'] = bridge_dur.sum() / duration if duration > 0 else np.nan

    # -----------------------------------------------------
    # KDE of detected t2
    # -----------------------------------------------------

    valid = t2.dropna().to_numpy()

    _, peak_pos, peak_height = kde_peaks(valid, prominence=p['kde_prominence'])

    feat['n_valid'] = len(valid)
    feat['n_peaks'] = len(peak_pos)

    film_peaks = peak_pos < p['h_liq']
    liquid_peaks = peak_pos >= p['h_liq']

    feat['peak_film_mm'] = peak_pos[film_peaks].min() if film_peaks.any() else np.nan
    feat['peak_film_density'] = (
        peak_height[film_peaks][np.argmin(peak_pos[film_peaks])]
        if film_peaks.any() else np.nan
    )
    feat['peak_liquid_density'] = (
        peak_height[liquid_peaks].max() if liquid_peaks.any() else 0.0
    )
    feat['liquid_to_film_peak_ratio'] = (
        feat['peak_liquid_density'] / feat['peak_film_density']
        if film_peaks.any() else np.nan
    )
    feat['bimodality_coefficient'] = bimodality_coefficient(valid)

    # -----------------------------------------------------
    # Film statistics
    # -----------------------------------------------------

    film = t2[states == 'film'].to_numpy()

    if len(film) >= 5:
        feat['film_median_mm'] = np.median(film)
        feat['film_q05_mm'] = np.quantile(film, 0.05)
        feat['film_q95_mm'] = np.quantile(film, 0.95)
        feat['film_cv'] = np.std(film) / np.mean(film)
    else:
        for key in ['film_median_mm', 'film_q05_mm', 'film_q95_mm', 'film_cv']:
            feat[key] = np.nan

    feat['film_wave_freq_hz'] = _film_wave_frequency(df, states)

    feat.update(_film_estimate_features(df, p['h_liq']))

    feat.update(_film_jump_features(df, p['h_liq'], p['jump_lag'], p['jump_min_change']))

    feat.update(_gap_features(states, time, duration))

    return feat


def _film_jump_features(df, h_liq, lag, min_change, min_pairs=30):
    """
    Temporal character of the detected film, from pairs of detected samples
    ``lag`` seconds apart (a fixed time lag keeps the 100 Hz and kHz runs
    comparable):

    * film_jump_frac : share of pairs whose thickness changes by more than
                       ``min_change`` -> abrupt changes (film / bubble passing
                       at a random height, dispersed) vs. smooth waves
    * film_lag_corr  : correlation between the two samples of a pair
                       (smooth waves stay correlated, bubbles do not)
    """

    time = df['time'].to_numpy(dtype=float)
    t2 = df['t2'].to_numpy(dtype=float).copy()
    t2[t2 >= h_liq] = np.nan

    if len(t2) < 2:
        return {'film_jump_frac': np.nan, 'film_lag_corr': np.nan}

    step = max(1, int(round(lag / np.nanmedian(np.diff(time)))))

    a, b = t2[:-step], t2[step:]
    valid = ~np.isnan(a) & ~np.isnan(b)

    if valid.sum() < min_pairs:
        return {'film_jump_frac': np.nan, 'film_lag_corr': np.nan}

    a, b = a[valid], b[valid]

    return {
        'film_jump_frac': float(np.mean(np.abs(b - a) > min_change)),
        'film_lag_corr': float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else np.nan,
    }


def film_estimate(df, h_liq=H_LIQ):
    """
    Bottom film thickness per sample [mm], including films thinner than the
    detection limit:

    * detected interface below h_liq      -> t2
    * t2 missing but t1 above the threshold -> (t1 - 0.985) n_water / n_bk7
      (film merged with the wall reading, same correction as the processing)
    * otherwise                           -> NaN (no film information)

    Liquid bridges (t2 >= h_liq) are excluded. Uses the unclipped ``t1_raw``
    column when present (processed runs), else ``t1`` (raw runs).
    """

    t2 = df['t2'].to_numpy(dtype=float)
    t1 = df['t1_raw' if 't1_raw' in df.columns else 't1'].to_numpy(dtype=float)

    merged = np.clip(np.nan_to_num(t1, nan=0) - T1_THRESHOLD, 0, None) * N_WATER / N_BK7

    estimate = np.where(
        ~np.isnan(t2),
        np.where(t2 < h_liq, t2, np.nan),
        np.where(merged > 0, merged, np.nan),
    )

    return estimate


def _film_estimate_features(df, h_liq):
    """
    Mean, standard deviation and coefficient of variation of the bottom film
    (``film_estimate``) and the fraction of samples where it is available.
    Mean and std are given in mm (``_mm``) and divided by CHANNEL_SIDE (``_rel``).
    """

    estimate = film_estimate(df, h_liq)
    valid = estimate[~np.isnan(estimate)]

    if len(valid) < 5:
        return {
            'film_est_coverage': len(valid) / len(estimate) if len(estimate) else 0.0,
            'film_est_mean_mm': np.nan,
            'film_est_std_mm': np.nan,
            'film_est_mean_rel': np.nan,
            'film_est_std_rel': np.nan,
            'film_est_cv': np.nan,
        }

    mean = valid.mean()
    std = valid.std()

    return {
        'film_est_coverage': len(valid) / len(estimate),
        'film_est_mean_mm': mean,
        'film_est_std_mm': std,
        'film_est_mean_rel': mean / CHANNEL_SIDE,
        'film_est_std_rel': std / CHANNEL_SIDE,
        'film_est_cv': std / mean if mean > 0 else np.nan,
    }


def relative_film_columns(df):
    """
    ``film_est_mean_rel`` and ``film_est_std_rel`` of a feature table; derived
    from the mm columns for tables of earlier versions that lack them.
    """

    if {'film_est_mean_rel', 'film_est_std_rel'} <= set(df.columns):
        return df['film_est_mean_rel'], df['film_est_std_rel']

    return df['film_est_mean_mm'] / CHANNEL_SIDE, df['film_est_std_mm'] / CHANNEL_SIDE


def _gap_features(states, time, duration, long_gap=GAP_LONG):
    """
    Structure of the missing t2 (undetected) samples:

    * gap_freq_hz        : number of undetected gaps per second
    * gap_mean_ms        : mean gap duration
    * gap_long_frac      : fraction of the undetected time in gaps >= long_gap
                           (long blackouts: film below detection for long periods)
    * gap_near_interface : fraction of gaps that touch a bubble/liquid sample
                           (missing at bubble noses/tails and slug fronts)
                           rather than being surrounded by wall film
    """

    state_values = states.to_numpy()
    segments, durations = _segments(state_values == 'undetected', time)

    if not segments:
        return {
            'gap_freq_hz': 0.0,
            'gap_mean_ms': 0.0,
            'gap_long_frac': 0.0,
            'gap_near_interface': 0.0,
        }

    tall_neighbour = []

    for start, end in segments:

        neighbours = []

        if start > 0:
            neighbours.append(state_values[start - 1])

        if end < len(state_values) - 1:
            neighbours.append(state_values[end + 1])

        tall_neighbour.append(any(s in ('bubble', 'liquid') for s in neighbours))

    return {
        'gap_freq_hz': len(segments) / duration if duration > 0 else np.nan,
        'gap_mean_ms': 1000 * durations.mean(),
        'gap_long_frac': durations[durations >= long_gap].sum() / durations.sum(),
        'gap_near_interface': float(np.mean(tall_neighbour)),
    }


def _film_wave_frequency(df, states):
    """
    Mean-crossing rate of the film thickness inside contiguous film
    segments, divided by two (one wave = two crossings) [Hz].
    """

    time = df['time'].to_numpy()
    t2 = df['t2'].to_numpy()

    segments, durations = _segments(states == 'film', time)

    if not segments:
        return np.nan

    film_mean = np.nanmean(t2[states.to_numpy() == 'film'])

    crossings = 0
    film_time = 0.0

    for (start, end), duration in zip(segments, durations):

        if end - start < 2:
            continue

        sign = np.sign(t2[start:end + 1] - film_mean)
        crossings += np.count_nonzero(np.diff(sign) != 0)
        film_time += duration

    return crossings / 2 / film_time if film_time > 0 else np.nan

# =========================================================
# Criterion
# =========================================================

def classify_frame(df, params=None):
    """
    Hierarchical objective flow pattern criterion, applied to a feature
    table (one row per run). ``params`` overrides entries of
    ``CRITERION_PARAMS``. Returns a numpy array of pattern names.

    1. Intermittent family: the liquid bridges the channel (t2 >= h_liq)
       at least f_bridge_hz times per second.
         - bubbly : interfaces seen are mostly bubbles, not a wall film
         - plug   : long (>= plug_min_duration) liquid bridges occupying
                    >= plug_min_f_liquid of the time
         - slug   : otherwise (short, aerated liquid slugs)
    2. Separated family: no liquid bridges. Based on the mean and standard
       deviation of the bottom film (``film_estimate``, including the film
       below the detection limit):
         - annular    : thin and uniform film (mean / H < annular_max_mean_rel and
                        std / H < annular_max_std_rel). Runs without any film
                        information also fall here (film always below
                        detection)
         - agitated film (churn / dispersed): liquid film present with
                        abruptly varying thickness (the probe sees either the
                        film or a bubble passing at a random height): mean / H >=
                        dispersed_min_mean_rel, film_jump_frac >= dispersed_min_jump
                        and the film detected in >= dispersed_min_detected of
                        the samples, or a KDE with >= dispersed_min_peaks peaks
                        (film + bubble layer); AND a mean pressure gradient at
                        least agitated_min_homog_ratio times the homogeneous
                        prediction (phases strongly mixed; below it the gas
                        slips over a separated layer and the run is wavy).
         - liquid-rich film (churn / dispersed): the bottom probe sees only a
                        thin, broken film (mean / H < liquid_rich_max_film_rel; the
                        interface is lost to aeration/scattering), but the
                        mean pressure gradient says the flow is mixed
                        (ratio_homog >= agitated_min_homog_ratio) and
                        liquid-dominated (phi2_l <= liquid_rich_max_phi2,
                        below every thin-film wavy/annular run). High-velocity
                        runs at the churn-annular transition.
           Agitated and liquid-rich runs are then split:
              - dispersed : KDE with >= dispersed_min_peaks peaks, or two-phase
                            multiplier phi2_l <= dispersed_max_phi2 (friction
                            dominated by the liquid)
              - churn     : otherwise (the gas adds much of the friction)
         - wavy       : otherwise (film with smooth waves; includes the
                        visual stratified runs, see ``REFERENCE_MAP``)

    Needs the pressure features (``src.pressure_features``) merged into
    ``df``: ratio_homog and phi2_l.
    """

    p = merge_params(CRITERION_PARAMS, params, 'CRITERION_PARAMS')

    bridge = df['bridge_freq_hz'].fillna(0).to_numpy() >= p['f_bridge_hz']
    bubble = df['bubble_share'].fillna(0).to_numpy() >= p['bubble_share_min']
    plug = (
        (df['bridge_mean_duration_s'].fillna(0).to_numpy() >= p['plug_min_duration'])
        & (df['f_liquid'].fillna(0).to_numpy() >= p['plug_min_f_liquid'])
    )

    # Bottom film mean and std divided by CHANNEL_SIDE (dimensionless)
    film_mean_rel, film_std_rel = relative_film_columns(df)
    film_mean = film_mean_rel.fillna(0).to_numpy()
    film_std = film_std_rel.fillna(0).to_numpy()

    annular = (film_mean < p['annular_max_mean_rel']) & (film_std < p['annular_max_std_rel'])

    multi_peak = df['n_peaks'].fillna(0).to_numpy() >= p['dispersed_min_peaks']
    jumpy_film = (
        (film_mean >= p['dispersed_min_mean_rel'])
        & (df['film_jump_frac'].fillna(0).to_numpy() >= p['dispersed_min_jump'])
        & (1 - df['f_undetected'].fillna(1).to_numpy() >= p['dispersed_min_detected'])
    )

    mixed = df['ratio_homog'].fillna(0).to_numpy() >= p['agitated_min_homog_ratio']
    phi2 = df['phi2_l'].fillna(np.inf).to_numpy()

    agitated = (multi_peak | jumpy_film) & mixed

    liquid_rich = (
        ~annular
        & ~agitated
        & mixed
        & (phi2 <= p['liquid_rich_max_phi2'])
        & (film_mean < p['liquid_rich_max_film_rel'])
    )
    agitated = agitated | liquid_rich

    liquid_friction = phi2 <= p['dispersed_max_phi2']
    dispersed = agitated & (multi_peak | liquid_friction)

    return np.select(
        [
            bridge & bubble,
            bridge & plug,
            bridge,
            annular,
            dispersed,
            agitated,
        ],
        ['bubbly', 'plug', 'slug', 'annular', 'dispersed', 'churn'],
        default='wavy',
    )


def classify_objective(feat, params=None):
    """
    Objective flow pattern of one run (``feat`` is a dict or Series of
    features). See ``classify_frame`` for the rules.
    """

    return classify_frame(pd.DataFrame([dict(feat)]), params)[0]

# =========================================================
# Tables
# =========================================================

def features_table(run_ids, base_path, trim=False, params=None, raw_path=None):
    """
    Feature table for a list of runs (one row per run). ``raw_path`` adds
    the unclipped t1 of the raw runs (needed for processed runs, see
    ``load_run``).
    """

    rows = []

    for run_id in run_ids:

        df = load_run(run_id, base_path, trim=trim, raw_path=raw_path)
        rows.append({'run_id': run_id, **run_features(df, params)})

    return pd.DataFrame(rows)


def balanced_accuracy(visual, objective, classes):
    """
    Mean recall over ``classes`` (runs whose visual label is not in
    ``classes`` are ignored). Returns (balanced accuracy, plain agreement).
    """

    visual = np.asarray(visual)
    objective = np.asarray(objective)

    mask = np.isin(visual, classes)
    hit = visual[mask] == objective[mask]

    recalls = [hit[visual[mask] == c].mean() for c in classes if (visual[mask] == c).any()]

    return float(np.mean(recalls)), float(hit.mean())
