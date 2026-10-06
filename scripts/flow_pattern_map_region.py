import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

import pandas as pd

import matplotlib.pyplot as plt

from src.plot_flow_pattern_maps import map_garimella_region
from src.flow_pattern_criteria import JL_MIN

silver_path = project_root / 'data/silver_layer/'
gold_path = project_root / 'data/gold_layer/'

evaluated_region_path = gold_path / 'flow_pattern_map_evaluated_region.pdf'

df_experimental = pd.read_pickle(silver_path / 'experimental_results.pkl')

plt.rcParams.update({
    "text.usetex": True,

    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],

    "ytick.color": "black",
    "xtick.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",

    "axes.linewidth": 0.8,

    "xtick.direction": "in",
    "ytick.direction": "in",
})

# The evaluated region covers the runs used in the objective flow pattern
# analysis: jl >= JL_MIN (the runs below it, mostly stratified, are excluded)
df_evaluated = df_experimental[df_experimental['jl'] >= JL_MIN]

print(f'{len(df_evaluated)} of {len(df_experimental)} runs (jl >= {JL_MIN} m/s)')

fig, ax = map_garimella_region(df_evaluated)

fig.savefig(evaluated_region_path, format='pdf', bbox_inches='tight')

print(f'Saved {evaluated_region_path}')
