"""
Independent checks of CLAUDE_PAPER_REVISION_BRIEF_V3.md on the reported results version.

Writes CSV files to data/gold_layer/<version>/verification_v3/ and prints a summary.
Run from the repository root:  python scripts/verify_v3.py [objective_criteria_10]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from src.dpdz_correlations import friction_churchill  # noqa: E402
from src.flow_pattern_criteria import (  # noqa: E402
    FEATURE_PARAMS, OBJECTIVE_PATTERNS, balanced_accuracy, classify_frame, features_table, load_run, sample_states,
)
from src.pressure_features import SQUARE_DUCT_PO, pressure_features, two_phase_state  # noqa: E402

version = sys.argv[1] if len(sys.argv) > 1 else 'objective_criteria_10'
gold = root / 'data/gold_layer' / version
out = gold / 'verification_v3'
out.mkdir(exist_ok=True)
silver = root / 'data/silver_layer'
RAW, V2 = silver / 'confocal_runs', silver / 'confocal_runs_processed_v2'

reported = pd.read_csv(gold / 'processed_v2/objective_flow_pattern.csv', index_col=0).set_index('run_id')
run_ids = reported.index.to_list()
exp = pd.read_pickle(silver / 'experimental_results.pkl').set_index('run_id').loc[run_ids].reset_index()


def section(title):
    print('\n' + '=' * 70 + '\n' + title + '\n' + '=' * 70)


# ---------------------------------------------------------------- 1. friction factor
section('1. Churchill-Jones friction factor (independent implementation)')


def churchill_darcy(re):
    """Churchill (1977), Darcy, smooth wall, written from the paper equations."""
    a = (2.457 * np.log(1 / ((7 / re) ** 0.9))) ** 16
    b = (37530 / re) ** 16
    return 8 * ((8 / re) ** 12 + (a + b) ** -1.5) ** (1 / 12)


re = np.array([100, 1000, 1500, 2000, 2299, 2300, 2301, 3000, 4000, 1e4, 1e5])
fric = pd.DataFrame({'Re': re, 'Re_star': 64 / 56.91 * re})
fric['f_independent'] = churchill_darcy(fric['Re_star'])
fric['f_src'] = friction_churchill(re, po=SQUARE_DUCT_PO)
fric['f_times_Re'] = fric['f_independent'] * re
fric['brief'] = fric['Re'].map({1000: 0.05691, 1500: 0.03794, 2000: 0.03027, 2300: 0.03721, 3000: 0.04267, 4000: 0.03914})
fric.to_csv(out / 'friction_table.csv', index=False)
print(fric.round(5).to_string(index=False))
assert np.allclose(fric['f_independent'], fric['f_src'], rtol=1e-12)
grid = np.linspace(500, 6000, 55001)
print('max relative step of f on a 0.1 Re grid, 500-6000:', f"{np.max(np.abs(np.diff(churchill_darcy(64 / 56.91 * grid))) / churchill_darcy(64 / 56.91 * grid)[1:]):.2e}")

# ---------------------------------------------------------------- 2. Reynolds windows and features
section('2. Reynolds numbers of the reference flows and pressure features')
state = two_phase_state(exp)
windows = []
for col in ('Re_l', 'Re_tp'):
    r = state[col]
    windows.append({'Re': col, 'below_2300': int((r < 2300).sum()), '2300_4000': int(((r >= 2300) & (r < 4000)).sum()),
                    'above_4000': int((r >= 4000).sum()), 'min': r.min(), 'max': r.max()})
windows = pd.DataFrame(windows)
windows.to_csv(out / 'reynolds_windows.csv', index=False)
print(windows.round(0).to_string(index=False))
pf = pressure_features(exp).set_index('run_id').loc[run_ids]
print('all R_h and phi_l^2 finite and positive:', bool(np.isfinite(pf[['ratio_homog', 'phi2_l']]).all().all() and (pf[['ratio_homog', 'phi2_l']] > 0).all().all()))
print('max |difference| to the reported features: R_h', float((pf['ratio_homog'] - reported['ratio_homog']).abs().max()),
      'phi_l^2', float((pf['phi2_l'] - reported['phi2_l']).abs().max()))

# ---------------------------------------------------------------- 3. gas-line convention and P_tp
section('3. Gas-line convention and test-section pressure')
A = 8e-3 ** 2
gas = pd.DataFrame({
    'run_id': exp['run_id'],
    'jg_over_Qg_A': exp['jg'] / (exp['Qg'] / 60000 / A),
    'rho_meter': exp['Pg'] / (287.0530 * (exp['Tg'] + 273.15)),
    'Ptp_minus_3p5_dp': exp['Ptp'] - 3.5 * exp['dp_exp'],
    'dpdz_over_dp_L': exp['dpdz_tp'] / (exp['dp_exp'] / 0.75),
})
gas['Gg_code_over_rho_jg'] = state['Gg'].to_numpy() / (gas['rho_meter'] * exp['jg']).to_numpy()
gas['G_minus_Gl_Gg'] = (state['g'] - state['Gl'] - state['Gg']).to_numpy()
gas.to_csv(out / 'gas_convention.csv', index=False)
print(gas.drop(columns='run_id').describe().loc[['min', 'max']].round(4).to_string())
print('Pg (absolute) range [kPa]:', exp['Pg'].min() / 1e3, '-', exp['Pg'].max() / 1e3)
tol = 150 / 0.75 / exp['dpdz_tp']
print(f'200 Pa/m relative to the measured gradient: median {100 * tol.median():.1f} %, range {100 * tol.min():.1f}-{100 * tol.max():.1f} %, '
      f'runs above 50 %: {int((tol > 0.5).sum())}')

# ---------------------------------------------------------------- 4. reproduction of the classification
section('4. Reproduction of features, classes and agreement')
feats = features_table(run_ids, V2, raw_path=RAW).set_index('run_id').loc[run_ids]
cols = [c for c in feats.columns if c in reported.columns and pd.api.types.is_numeric_dtype(feats[c])]
diff = (feats[cols].astype(float) - reported[cols].astype(float)).abs().max()
print('features with |difference| > 1e-9:', diff[diff > 1e-9].to_dict())
table = feats.join(pf[['ratio_homog', 'phi2_l']])
pred = classify_frame(table.reset_index())
visual = reported['reduced_pattern'].to_numpy()
print('classes equal to the reported ones:', bool((pred == reported['objective_pattern'].to_numpy()).all()))
ba, ag = balanced_accuracy(visual, pred, OBJECTIVE_PATTERNS)
order = ['bubbly', 'plug', 'slug', 'dispersed', 'churn', 'wavy', 'annular']
cm = pd.crosstab(pd.Categorical(visual, order), pd.Categorical(pred, order), rownames=['visual'], colnames=['objective'], dropna=False)
cm.to_csv(out / 'confusion_matrix.csv')
print(cm.to_string())
print(f'agreement {100 * ag:.1f} % ({int((pred == visual).sum())}/{len(visual)}), balanced accuracy {100 * ba:.1f} %, '
      f'disagreements {int((pred != visual).sum())}')
pr = pd.DataFrame({c: {'recall': f'{int(((visual == c) & (pred == c)).sum())}/{int((visual == c).sum())}',
                       'precision': f'{int(((visual == c) & (pred == c)).sum())}/{int((pred == c).sum())}'} for c in order}).T
print(pr.to_string())

# ---------------------------------------------------------------- 5. optical inversion: origin of delta_f values
section('5. Origin of the bottom-film values (Eq. 12)')
rows = []
for run_id in run_ids:
    proc = load_run(run_id, V2, raw_path=RAW)
    raw = load_run(run_id, RAW, trim=True)
    t2p = proc['t2'].to_numpy(float)
    t2r = raw['t2'].to_numpy(float)
    keep = ~proc['mixer_wave'].astype(bool).to_numpy() if 'mixer_wave' in proc else np.ones(len(proc), bool)
    first = keep & ~np.isnan(t2p) & (t2p < FEATURE_PARAMS['h_liq'])
    t1 = proc['t1_raw'].to_numpy(float)
    second = keep & np.isnan(t2p) & (np.nan_to_num(t1) > 0.985)
    rows.append({
        'run_id': run_id, 'pattern': reported.loc[run_id, 'reduced_pattern'],
        'detected_no_wall_excess': int((first & ~np.isnan(t2r) & np.isclose(t2p, t2r, atol=1e-6)).sum()),
        'detected_plus_wall_excess': int((first & ~np.isnan(t2r) & ~np.isclose(t2p, t2r, atol=1e-6)).sum()),
        'substituted_excess': int((first & np.isnan(t2r)).sum()),
        'indirect_second_branch': int(second.sum()),
    })
origin = pd.DataFrame(rows)
origin.to_csv(out / 'film_value_origin.csv', index=False)
kinds = ['detected_no_wall_excess', 'detected_plus_wall_excess', 'substituted_excess', 'indirect_second_branch']
by = origin.groupby('pattern')[kinds].sum()
by.loc['all'] = by.sum()
print((100 * by.div(by.sum(axis=1), axis=0)).round(1).to_string())

# ---------------------------------------------------------------- 6. sampling: jump lag, bridges, mixer gaps
section('6. Sampling: jump lag per rate, one-sample bridges, gaps left by removed mixer waves')
lag_rows, gap_rows = [], []
for run_id in run_ids:
    proc = load_run(run_id, V2, raw_path=RAW)
    dt = float(np.nanmedian(np.diff(proc['time'])))
    step = max(1, int(round(FEATURE_PARAMS['jump_lag'] / dt)))
    lag_rows.append({'run_id': run_id, 'rate_khz': round(1 / dt / 1000, 2), 'step': step, 'lag_ms': 1000 * step * dt})
    if 'mixer_wave' in proc and proc['mixer_wave'].any():
        kept = proc.loc[~proc['mixer_wave'].astype(bool)].reset_index(drop=True)
        t = kept['time'].to_numpy()
        gap = np.diff(t) > 1.5 * dt
        liquid = (sample_states(kept) == 'liquid').to_numpy()
        bridge_across = int((gap & liquid[:-1] & liquid[1:]).sum())
        # jump pairs (step apart in index) whose time difference is not step * dt
        pair_dt = t[step:] - t[:-step]
        gap_rows.append({'run_id': run_id, 'removed': int(proc['mixer_wave'].sum()), 'gaps': int(gap.sum()),
                         'bridges_across_gap': bridge_across,
                         'jump_pairs_longer_lag_pct': 100 * float(np.mean(pair_dt > (step + 0.5) * dt))})
lags = pd.DataFrame(lag_rows)
lags.to_csv(out / 'jump_lag.csv', index=False)
print(lags.groupby('rate_khz').agg(runs=('run_id', 'size'), step=('step', 'first'), lag_ms=('lag_ms', 'first')).round(2).to_string())
gaps = pd.DataFrame(gap_rows)
gaps.to_csv(out / 'mixer_gaps.csv', index=False)
print(f"runs with removed mixer waves: {len(gaps)}; bridges joined across a removed gap: {int(gaps['bridges_across_gap'].sum())}; "
      f"jump pairs spanning a gap: max {gaps['jump_pairs_longer_lag_pct'].max():.2f} % of a run's pairs")
one = reported[(reported['bridge_count'] > 0) & (reported['sampling_hz'] < 200)]
print(f"0.1 kHz runs with bridges: {len(one)}; minimum mean bridge duration {1000 * one['bridge_mean_duration_s'].min():.1f} ms "
      f"(one sample = {1000 / one['sampling_hz'].median():.1f} ms)")
rec = reported[['bridge_count', 'duration_s', 'bridge_freq_hz']]
single = rec[(rec['bridge_count'] == 1) & (rec['bridge_freq_hz'] >= 0.1)]
print(f'runs above f_bridge = 0.1 Hz with a single bridge: {len(single)} (records {single.duration_s.min():.1f}-{single.duration_s.max():.1f} s)')

# ---------------------------------------------------------------- 7. KDE normalization (Fig. 19)
section('7. Integral of the Fig. 19 densities over the plotted range 0-9 mm')
KDE_RUNS = {'bubbly': [65, 61, 60, 59, 50, 53], 'plug': [35, 33, 36, 31, 28, 66], 'slug': [15, 14, 13, 25, 87, 56],
            'dispersed': [81, 80, 74, 75, 76, 77], 'churn': [26, 27, 94, 138, 1, 140], 'wavy': [45, 41, 148, 47, 46, 44]}
kde_rows = []
for pattern, ids in KDE_RUNS.items():
    for run_id in ids:
        t2 = load_run(run_id, V2)['t2'].dropna().to_numpy()
        k = gaussian_kde(t2, bw_method='scott')
        kde_rows.append({'pattern': pattern, 'run_id': run_id, 'n': len(t2), 'integral_0_9': k.integrate_box_1d(0, 9),
                         'mass_below_0': k.integrate_box_1d(-np.inf, 0)})
kde = pd.DataFrame(kde_rows)
kde.to_csv(out / 'kde_integrals.csv', index=False)
print(kde.groupby('pattern')[['integral_0_9', 'mass_below_0']].agg(['min', 'max']).round(4).to_string())
print('\nwritten to', out)
