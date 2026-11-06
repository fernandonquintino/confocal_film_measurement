"""
Evaluation of a convective scaling of the time-based rules of the objective
criterion (reproducibility of the Stage 1 analysis of the Fluids paper).

The reported criterion keeps its dimensional time thresholds (f_bridge >= 0.1 Hz,
t_bridge >= 20 ms). This script only evaluates the alternatives, with every other
rule and threshold unchanged:

  family split : St = f_bridge D_h / U       instead of f_bridge
  plug split   : L* = t_bridge U / D_h       instead of t_bridge

for U = j_m = jl + jg, U = jl, and U = U_t = C0 j_m + U_d (C0 = 1.2,
U_d = 0.54 sqrt(g D_h); Bendiksen-type assumption, not validated for square ducts).
Superficial velocities at the meter conditions. Thresholds are swept over the
midpoints of the observed values; scores are the balanced accuracy (BA) and the
agreement (AG) of the whole criterion against the visual classes, in-sample and
with the threshold chosen by leave-one-out (LOO).

Note: L* is NOT a liquid-slug length. A liquid bridge is the unaerated stretch of
a slug; entrained bubbles split aerated slugs into several short bridges.

Usage:  python scripts/convective_scaling_eval.py [objective_criteria_N]
Output: data/gold_layer/<version>/convective_scaling/*.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import src.flow_pattern_criteria as fc  # noqa: E402

G = 9.81
D_H = fc.CHANNEL_SIDE * 1e-3        # hydraulic diameter of the square channel [m]


def velocity_scales(df):
    jm = (df['jl'] + df['jg']).to_numpy()
    return {'j_m': jm, 'j_l': df['jl'].to_numpy(), 'U_t': 1.2 * jm + 0.54 * np.sqrt(G * D_H)}


def candidates(x):
    """Geometric midpoints of the positive observed values (the variables span decades)."""
    u = np.unique(x[np.isfinite(x) & (x > 0)])
    return np.concatenate([[u[0] * 0.5], np.sqrt(u[:-1] * u[1:]), [u[-1] * 1.01]])


def classify_with(df, column, values, key, threshold):
    """v10/v11 criterion with one rule input replaced by ``values``."""
    d = df.copy()
    d[column] = values
    return fc.classify_frame(d, {key: threshold})


RULES = {   # rule name -> (feature column replaced, threshold key)
    'family': ('bridge_freq_hz', 'f_bridge_hz'),
    'plug': ('bridge_mean_duration_s', 'plug_min_duration'),
}


def evaluate(df, visual, rule, name, values):
    column, key = RULES[rule]
    base = df['objective_pattern'].to_numpy()
    thresholds = candidates(values)
    preds = np.array([classify_with(df, column, values, key, t) for t in thresholds])

    sweep = pd.DataFrame({
        'rule': rule, 'variable': name, 'threshold': thresholds,
        'balanced_accuracy': [fc.balanced_accuracy(visual, p, fc.OBJECTIVE_PATTERNS)[0] for p in preds],
        'agreement': [fc.balanced_accuracy(visual, p, fc.OBJECTIVE_PATTERNS)[1] for p in preds],
        'changed_vs_reported': [(p != base).sum() for p in preds],
    })

    # Leave-one-out: threshold (middle of the best plateau) chosen without the left-out run
    loo = np.empty(len(df), dtype=object)
    for i in range(len(df)):
        keep = np.arange(len(df)) != i
        scores = np.array([fc.balanced_accuracy(visual[keep], p[keep], fc.OBJECTIVE_PATTERNS)[0] for p in preds])
        best = np.flatnonzero(np.isclose(scores, scores.max()))
        loo[i] = preds[best[len(best) // 2], i]

    bmax = sweep['balanced_accuracy'].max()
    plateau = sweep[np.isclose(sweep['balanced_accuracy'], bmax)]
    near = sweep[sweep['balanced_accuracy'] >= bmax - 0.01]
    best = plateau.iloc[len(plateau) // 2]
    loo_ba, loo_ag = fc.balanced_accuracy(visual, loo, fc.OBJECTIVE_PATTERNS)

    summary = {
        'rule': rule, 'variable': name,
        'best_threshold': best['threshold'], 'best_balanced_accuracy': bmax, 'best_agreement': best['agreement'],
        'near_optimal_min': near['threshold'].min(), 'near_optimal_max': near['threshold'].max(),
        'near_optimal_ratio': near['threshold'].max() / near['threshold'].min(),
        'changed_vs_reported': int(best['changed_vs_reported']),
        'loo_balanced_accuracy': loo_ba, 'loo_agreement': loo_ag,
    }
    return sweep, summary


def main(version):
    gold = ROOT / 'data/gold_layer' / version
    df = pd.read_pickle(gold / 'processed_v2/objective_flow_pattern.pkl')
    visual = df['reduced_pattern'].to_numpy()
    out = gold / 'convective_scaling'
    out.mkdir(exist_ok=False)

    f = df['bridge_freq_hz'].fillna(0).to_numpy()
    t = df['bridge_mean_duration_s'].fillna(0).to_numpy()
    U = velocity_scales(df)

    variables = [('family', 'f_bridge [Hz] (reported)', f), ('plug', 't_bridge [s] (reported)', t)]
    variables += [('family', f'St = f_bridge D_h/{k}', f * D_H / u) for k, u in U.items()]
    variables += [('plug', f'L* = t_bridge {k}/D_h', t * u / D_H) for k, u in U.items()]

    sweeps, summaries = zip(*(evaluate(df, visual, *v) for v in variables))
    pd.concat(sweeps).to_csv(out / 'threshold_sweeps.csv', index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / 'summary.csv', index=False)

    # Per-run values, and the runs changed when both rules use j_m (in-sample best thresholds)
    per_run = pd.DataFrame({
        'run_id': df['run_id'], 'visual': visual, 'reported': df['objective_pattern'],
        'jl': df['jl'], 'jg': df['jg'], 'sampling_hz': df['sampling_hz'],
        'f_bridge_hz': f, 't_bridge_s': t,
        **{f'St_{k}': f * D_H / u for k, u in U.items()},
        **{f'Lstar_{k}': t * u / D_H for k, u in U.items()},
    })
    s = summary.set_index('variable')['best_threshold']
    d = df.copy()
    d['bridge_freq_hz'] = per_run['St_j_m'].to_numpy()
    d['bridge_mean_duration_s'] = per_run['Lstar_j_m'].to_numpy()
    per_run['both_convective_j_m'] = fc.classify_frame(
        d, {'f_bridge_hz': s['St = f_bridge D_h/j_m'], 'plug_min_duration': s['L* = t_bridge j_m/D_h']})
    per_run.to_csv(out / 'per_run.csv', index=False)

    # Signal-processing time scales in convective form: shorter than one sample?
    dt = 1 / df['sampling_hz'].to_numpy()
    slow = dt > 0.005
    rows = []
    for name, tau in [('min_dwell', fc.FEATURE_PARAMS['min_dwell']), ('jump_lag', fc.FEATURE_PARAMS['jump_lag'])]:
        c = tau * np.median(U['j_m']) / D_H
        tau_run = c * D_H / U['j_m']
        rows.append({'parameter': name, 'value_s': tau, 'dimensionless_at_median_j_m': c,
                     'run_value_min_s': tau_run.min(), 'run_value_max_s': tau_run.max(),
                     'runs_below_one_sample': int((tau_run < dt).sum()),
                     'runs_100hz_below_one_sample': int(((tau_run < dt) & slow).sum()), 'runs_100hz': int(slow.sum())})
    pd.DataFrame(rows).to_csv(out / 'processing_time_scales.csv', index=False)

    print(summary.round(4).to_string(index=False))
    print('Output:', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'objective_criteria_11')
