import pandas as pd

from src.combined_properties import properties_df
from src.dimensionless_number import Reynolds
from src.dpdz_correlations import dpdz_mcaddams_1942, friction_f

# =========================================================
# Mean pressure drop features
# =========================================================
#
# Only the mean pressure gradient of each run is available (no dynamic
# pressure signal). It is normalised by two reference gradients:
#
# * ratio_homog : measured / homogeneous (McAdams, no slip) prediction.
#                 >= 1 when the phases are strongly mixed (churn, dispersed);
#                 < 1 when the gas slips over a separated liquid layer.
# * phi2_l      : measured / liquid-alone gradient (Lockhart-Martinelli
#                 two-phase multiplier). Small when the friction is
#                 dominated by the liquid (dispersed), large when the gas
#                 adds much of it (churn, separated flows).

PRESSURE_FEATURES = ['ratio_homog', 'phi2_l']


def pressure_features(df_experimental, d=8e-3):
    """
    Pressure-drop features per run from the experimental table (needs
    run_id, Ttp, Ptp, jl, jg and dpdz_tp).
    """

    props = properties_df(df_experimental, d=d)

    rows = []

    for _, r in props.iterrows():

        dpdz_homog = dpdz_mcaddams_1942(r['g'], r['d'], r['x'], r['rhol'], r['rhog'], r['mul'], r['mug'])

        re_l = Reynolds(r['Gl'], r['d'], r['mul'])
        dpdz_l = r['Gl'] ** 2 / (2 * r['d'] * r['rhol']) * friction_f(re_l)

        rows.append({'run_id': r['run_id'], 'dpdz_homog': dpdz_homog, 'dpdz_l': dpdz_l})

    out = pd.DataFrame(rows).merge(df_experimental[['run_id', 'dpdz_tp']], on='run_id')

    out['ratio_homog'] = out['dpdz_tp'] / out['dpdz_homog']
    out['phi2_l'] = out['dpdz_tp'] / out['dpdz_l']

    return out[['run_id', 'dpdz_homog', 'dpdz_l'] + PRESSURE_FEATURES]
