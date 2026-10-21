"""
Final reconciliation of the manuscript numbers with the current results (CONFOCAL_PAPER_V6_FINAL_REVIEW_CHECKLIST.md, P1-D).

Complements scripts/verify_v3.py (features, classes, confusion matrix, friction factor, optical shares), which
should be run first. Prints PASS/FAIL per check. Run from the repository root:
    python scripts/verify_final_v6.py [objective_criteria_10]
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from src.flow_pattern_criteria import CRITERION_PARAMS  # noqa: E402

version = sys.argv[1] if len(sys.argv) > 1 else 'objective_criteria_10'
gold = root / 'data/gold_layer' / version
art = root / 'artigo'
tex = (art / 'confocal_paper.tex').read_text(encoding='utf-8')
runs = pd.read_csv(gold / 'processed_v2/objective_flow_pattern.csv', index_col=0)
order = ['bubbly', 'plug', 'slug', 'dispersed', 'churn', 'wavy', 'annular']
results = []


def check(name, ok, detail=''):
    results.append((name, 'PASS' if ok else 'FAIL', detail))


def table_rows(path):
    """Rows of a LaTeX table body as lists of cells (between midrule and bottomrule)."""
    body = (art / path).read_text(encoding='utf-8').split('\\midrule', 1)[1].split('\\bottomrule')[0]
    return [[c.strip() for c in r.split('&')] for r in body.split('\\\\') if '&' in r]


# 1. Table 9 medians ------------------------------------------------------------
features = [('bridge_freq_hz', 1), ('bridge_mean_duration_s', 1000), ('f_liquid', 1), ('bubble_share', 1), ('n_peaks', 1),
            ('film_est_mean_mm', 1), ('film_est_std_mm', 1), ('film_jump_frac', 1), ('f_undetected', 1),
            ('ratio_homog', 1), ('phi2_l', 1)]
bad = []
for (col, scale), row in zip(features, table_rows('tables/feature_medians.tex')):
    for pattern, cell in zip(order, row[1:]):
        median = scale * runs.loc[runs['reduced_pattern'] == pattern, col].median()
        if cell == '--':
            ok = median == 0 or (col == 'n_peaks' and pattern == 'annular')
        else:
            decimals = len(cell.split('.')[1]) if '.' in cell else 0
            ok = abs(float(cell) - median) <= 0.5 * 10 ** -decimals + 1e-9
        if not ok:
            bad.append(f'{col}/{pattern}: table {cell}, data {median:.4g}')
check('Table 9 medians (77 cells)', not bad, '; '.join(bad))

# 2. thresholds -----------------------------------------------------------------
manifest = json.loads((gold / 'run_manifest.json').read_text())
check('Run manifest thresholds = CRITERION_PARAMS', manifest['criterion_params'] == CRITERION_PARAMS)
p = CRITERION_PARAMS
expected = {
    'f_bridge 0.1 Hz': p['f_bridge_hz'] == 0.1, 'S_b 0.5': p['bubble_share_min'] == 0.5,
    't_bridge 20 ms': p['plug_min_duration'] == 0.02, 'f_liquid 0.10': p['plug_min_f_liquid'] == 0.1,
    'annular 0.1/0.2 mm': (p['annular_max_mean'], p['annular_max_std']) == (0.1, 0.2),
    'R_h 1': p['agitated_min_homog_ratio'] == 1.0, 'dispersed phi 5.5': p['dispersed_max_phi2'] == 5.5,
    'liquid-rich phi 20': p['liquid_rich_max_phi2'] == 20.0, 'liquid-rich film 0.5 mm': p['liquid_rich_max_film'] == 0.5,
    'agitated N_p 2, mean 0.4, J 0.4, detected 15 %': (p['dispersed_min_peaks'], p['dispersed_min_mean'], p['dispersed_min_jump'],
                                                       p['dispersed_min_detected']) == (2, 0.4, 0.4, 0.15),
}
check('Code thresholds = values stated in the paper', all(expected.values()), ', '.join(k for k, v in expected.items() if not v))
tree = (art / 'figures/paper_figures/fig11_decision_tree.tex').read_text(encoding='utf-8')
rules = (art / 'tables/criterion_rules.tex').read_text(encoding='utf-8')
need = {
    'Fig. 11': (tree, [r'f_{bridge} \geq 0.1', r'S_b \geq 0.5', r'\bar t_{bridge} \geq 20', r'f_{liquid} \geq 0.10',
                       r'\bar\delta_f < 0.1', r'\sigma_f < 0.2', r'R_h \geq 1', r'\phi_l^2 \leq 5.5', r'\phi_l^2 \leq 20',
                       r'\bar\delta_f < 0.5', r'J \geq 0.4']),
    'Table 8': (rules, [r'f_{bridge} \geq 0.1', r'S_b \geq 0.5', r'\bar t_{bridge} \geq 20', r'f_{liquid} \geq 0.10',
                        r'\bar\delta_f < 0.1', r'\sigma_f < 0.2', r'R_h \geq 1', r'\phi_l^2 \leq 5.5', r'\phi_l^2 \leq 20',
                        r'\bar\delta_f < 0.5', r'J \geq 0.4']),
    '§3.4 text': (tex, [r'$S_b \geq 0.5$', r'\bar t_{bridge} \geq \SI{20}{\milli\second}', r'f_{liquid} \geq 0.10',
                        r'$R_h \geq 1$', r'\phi_l^2 \leq 5.5', r'\phi_l^2 \leq 20', r'$J \geq 0.4$']),
}
for where, (text, items) in need.items():
    flat = re.sub(r'\s+', ' ', text.replace('~', ' ').replace('\\,', ' ').replace('$', '$'))
    missing = [i for i in items if re.sub(r'\s+', ' ', i) not in flat and i.replace(' ', '') not in flat.replace(' ', '')]
    check(f'Thresholds printed in {where}', not missing, ', '.join(missing))
sweep = pd.read_csv(gold / 'processed_v2/sensitivity_thresholds.csv')
best = {k: sweep.loc[(sweep.parameter == k) & np.isclose(sweep.value, p[k]), 'balanced_accuracy'] for k in
        ['f_bridge_hz', 'bubble_share_min', 'plug_min_duration', 'annular_max_mean', 'agitated_min_homog_ratio',
         'dispersed_max_phi2', 'liquid_rich_max_phi2', 'liquid_rich_max_film']}
check('Fig. 24: every chosen value is on its sweep and gives BA 88.1 %',
      all(len(v) == 1 and round(100 * v.iloc[0], 1) == 88.1 for v in best.values()),
      ', '.join(f'{k}: {list((100 * v).round(1))}' for k, v in best.items() if not (len(v) == 1 and round(100 * v.iloc[0], 1) == 88.1)))
nb = (root / 'notebooks/paper_figures.ipynb').read_text(encoding='utf-8')
check('paper_figures.ipynb reads the same results version', f"CRITERIA_VERSION = '{version}'" in nb)

# 3. Tables 10 and 11 ----------------------------------------------------------
var = pd.read_csv(gold / 'variant_comparison_criterion.csv').set_index('data')
t10 = {r[0]: (float(r[1].split('\\')[0]), float(r[2].split('\\')[0])) for r in table_rows('tables/variant_comparison.tex')}
pairs = {'Raw': 'regular', 'Raw, trimmed only': 'regular (trimmed only)',
         'Processed (trimming and tube-thickness correction)': 'processed', 'Processed and mixer waves removed': 'processed_v2'}
bad = [k for k, v in pairs.items() if t10[k] != (round(100 * var.loc[v, 'agreement'], 1), round(100 * var.loc[v, 'balanced_accuracy'], 1))]
check('Table 10 = variant_comparison_criterion.csv', not bad, ', '.join(bad))

dis = runs[runs.reduced_pattern != runs.objective_pattern]
bad = []
counts = {}
for r in table_rows('tables/disagreements.tex'):
    vis, obj = [s.strip().lower() for s in r[0].split('$\\rightarrow$')]
    counts.setdefault((vis, obj), []).append((int(r[1]), r[2], r[3]))
for (vis, obj), rows in counts.items():
    g = dis[(dis.reduced_pattern == vis) & (dis.objective_pattern == obj)]
    if sum(n for n, _, _ in rows) != len(g):
        bad.append(f'{vis}->{obj}: table {sum(n for n, _, _ in rows)}, data {len(g)}')
    if len(rows) == 1:
        n, jl, jg = rows[0]
        rng = lambda s: f'{s.min():.2f}--{s.max():.2f}' if round(s.min(), 2) != round(s.max(), 2) else f'{s.min():.2f}'
        if (jl, jg) != (rng(g.jl), rng(g.jg)):
            bad.append(f'{vis}->{obj}: ranges table {jl}, {jg}; data {rng(g.jl)}, {rng(g.jg)}')
check('Table 11: 19 runs, counts and velocity ranges per group',
      len(dis) == 19 and sum(n for v in counts.values() for n, _, _ in v) == 19 and not bad, '; '.join(bad))

# 4. numbers in the text --------------------------------------------------------
agree = (runs.reduced_pattern == runs.objective_pattern).mean()
check('Abstract / §4.5 / §5: 87.9 % (138 of 157), 88.1 %, raw 76.4 %',
      all(s in tex for s in [r'\SI{87.9}{\percent}', '138 of 157', r'\SI{88.1}{\percent}', r'\SI{76.4}{\percent}'])
      and round(100 * agree, 1) == 87.9)
ph = pd.read_csv(gold / 'verification/downsample_100hz_phases.csv')
per = ph.assign(c=ph.objective_100hz != ph.objective_khz).groupby('offset').c.sum()
check('§4.6 downsampling: 1 to 13 of 59 runs', ph.run_id.nunique() == 59 and (per.min(), per.max()) == (1, 13)
      and '1 to 13 of them change class' in tex)
pb = pd.read_csv(gold / 'verification/pressure_bounds.csv')
n_dp, n_all = int(pb.patterns_dp_only.str.contains('/').sum()), int(pb.patterns_all.str.contains('/').sum())
check('§4.6 pressure tolerance: 2 runs (Δp), 7 runs (worst case)', (n_dp, n_all) == (2, 7)
      and 'two runs could change' in tex and 'and seven,' in tex, f'data {n_dp}, {n_all}')

# 5. §3.2 shares ----------------------------------------------------------------
o = pd.read_csv(gold / 'verification_v3/film_value_origin.csv')
k = ['detected_no_wall_excess', 'detected_plus_wall_excess', 'substituted_excess', 'indirect_second_branch']
tot = o[k].sum()
first = tot[k[:3]].sum()
ann = o[o.pattern == 'annular'][k].sum()
shares = [100 * tot[k[0]] / first, 100 * tot[k[1]] / first, 100 * tot[k[2]] / first,
          100 * tot[k[3]] / tot.sum(), 100 * ann[k[3]] / ann.sum()]
check('§3.2 shares 90 / 9 / 0.7 % of first-branch values; 39 % of all values; 98 % in annular',
      [round(shares[0]), round(shares[1]), round(shares[2], 1), round(shares[3]), round(shares[4])] == [90, 9, 0.7, 39, 98],
      ', '.join(f'{s:.1f}' for s in shares))

for name, status, detail in results:
    print(f'{status}  {name}' + (f'  [{detail}]' if detail else ''))
print(f"\n{sum(s == 'PASS' for _, s, _ in results)}/{len(results)} checks passed")
