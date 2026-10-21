"""Minimum measured film thickness (Table 4 of the paper) and the share of missing second thickness
that the tube-thickness correction would fill, on the trimmed raw runs of the 157 analysed runs.

Same definitions as notebooks/thickness.ipynb (calculate_t2_stats and the threshold summary), but:
* each raw run is trimmed to start at the first valid first thickness (as in the data treatment),
  so the missing rate no longer counts the initial block without signal;
* only the 157 runs with jl >= JL_MIN;
* refractive indexes of the treatment (n_water = 1.333, n_BK7 = 1.519).

Run from the project root:  python scripts/minimum_film_thickness_table.py
Writes artigo/tables/minimum_film_thickness.tex and data/gold_layer/minimum_film_thickness.csv,
and prints the fill percentages quoted in the text.
"""
import sys
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.flow_pattern_criteria import JL_MIN, N_BK7, N_WATER, T1_THRESHOLD, load_run, reference_pattern  # noqa: E402

MIN_T2 = 0.373
N_ROWS = 5

silver = project_root / 'data/silver_layer'
df_exp = pd.read_pickle(silver / 'experimental_results.pkl')
df_exp = df_exp[df_exp['jl'] >= JL_MIN].copy()
df_exp['reference'] = reference_pattern(df_exp['reduced_pattern'])

rows = []
fill = {thr: [0, 0] for thr in (0.0, 0.1, 0.2, 0.3, MIN_T2)}

for run_id, pattern in zip(df_exp['run_id'], df_exp['reference']):
    df = load_run(run_id, silver / 'confocal_runs', trim=True)
    t2 = df['t2'].dropna()
    rows.append({
        'run_id': run_id, 'pattern': pattern,
        'missing_rate_t2': 100 * df['t2'].isna().mean(),
        'min_t2': t2.min(), 'q05_t2': t2.quantile(0.05), 'q10_t2': t2.quantile(0.10),
    })
    d2 = (df['t1'] - T1_THRESHOLD).clip(lower=0) * N_WATER / N_BK7
    missing = df['t2'].isna()
    for thr in fill:
        fill[thr][0] += int(missing.sum())
        fill[thr][1] += int((missing & (d2 > thr)).sum())

stats = pd.DataFrame(rows)
stats.to_csv(project_root / 'data/gold_layer/minimum_film_thickness.csv', index=False)

stable = stats[stats['missing_rate_t2'] < 50].nsmallest(N_ROWS, 'q05_t2')
unstable = stats[stats['missing_rate_t2'] >= 50].nsmallest(N_ROWS, 'q05_t2')
print('stable runs (< 50 % missing):\n', stable.round(3).to_string(index=False))
print('unstable runs (>= 50 % missing):\n', unstable.round(3).to_string(index=False))
print('runs with < 50 % missing:', int((stats['missing_rate_t2'] < 50).sum()), 'of', len(stats))
for thr, (n_missing, n_filled) in fill.items():
    print(f'threshold {thr:.3f} mm: fills {100 * n_filled / n_missing:.2f} % of the missing second thickness')


def table_rows(df):
    return '\n'.join(f'{r.pattern.capitalize():<10}& {r.q05_t2:.3f} & {r.q10_t2:.3f} & {r.min_t2:.3f} & {r.missing_rate_t2:.2f} \\\\'
                     for r in df.itertuples())


tex = r'''\begin{table}[!htbp]
\centering
\small
\caption{Minimum and low-percentile measured film thicknesses for selected runs - values in \si{\milli\meter}}
\label{tab:minimum_thickness}
\begin{tabular}{lAAAA}
\toprule
Flow pattern & q_{05} & q_{10} & \text{Min} & \text{Missing rate (\%)} \\
\midrule
''' + table_rows(stable) + '\n\\midrule\n' + table_rows(unstable) + r'''
\bottomrule
\end{tabular}
\end{table}
'''
(project_root / 'artigo/tables/minimum_film_thickness.tex').write_text(tex, encoding='utf-8')
print('written artigo/tables/minimum_film_thickness.tex')
