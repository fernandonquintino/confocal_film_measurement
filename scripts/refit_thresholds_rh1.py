"""Refit the thresholds on objective_criteria_9 features with the R_h limit fixed at 1 (coarse grid of kde_plots cell 41)."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from src.flow_pattern_criteria import CRITERION_PARAMS, OBJECTIVE_PATTERNS, balanced_accuracy, classify_frame  # noqa: E402

d = pd.read_csv(root / 'data/gold_layer/objective_criteria_9/processed_v2/objective_flow_pattern.csv', index_col=0)
grid = {
    'f_bridge_hz': [0.05, 0.1, 0.2, 0.5],
    'plug_min_duration': [0.01, 0.015, 0.02, 0.025, 0.03],
    'annular_max_mean': [0.05, 0.1, 0.15],
    'annular_max_std': [0.1, 0.2, 0.3],
    'dispersed_min_mean': [0.3, 0.4, 0.5, 0.6],
    'dispersed_min_jump': [0.3, 0.35, 0.4, 0.45],
    'dispersed_max_phi2': [4.0, 5.0, 5.5, 6.0],
    'liquid_rich_max_phi2': [12.0, 15.0, 18.0, 20.0],
}
keys = list(grid)
rows = []
for combo in itertools.product(*grid.values()):
    params = dict(zip(keys, combo))
    pred = classify_frame(d, params)
    ba, ag = balanced_accuracy(d['reduced_pattern'], pred, OBJECTIVE_PATTERNS)
    changed = sum(v != CRITERION_PARAMS[k] for k, v in params.items())
    rows.append({**params, 'ba': ba, 'ag': ag, 'n_changed': changed})
df = pd.DataFrame(rows)
best = df['ba'].max()
top = df[df['ba'] >= best - 1e-12].sort_values(['n_changed'])
print('combinations', len(df), 'best BA', round(100 * best, 2), 'ties', len(top))
print(top.head(8).to_string(index=False))
print('changed-parameter frequency among ties:')
for k in keys:
    print(' ', k, top[k].value_counts().to_dict())
