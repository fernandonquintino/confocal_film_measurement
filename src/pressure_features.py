import numpy as np
import pandas as pd

from src.combined_properties import properties_df
from src.dpdz_correlations import friction_churchill, friction_f

# =========================================================
# Mean pressure drop features
# =========================================================
#
# Only the mean pressure gradient of each run is available (no dynamic
# pressure signal). It is normalised by two frictional reference gradients:
#
# * ratio_homog : measured / homogeneous (McAdams, no slip) prediction.
#                 Observed >= 1 for churn and dispersed runs and < 1 for
#                 wavy runs.
# * phi2_l      : measured / liquid-alone gradient (liquid-reference
#                 two-phase multiplier). Small when the measured gradient
#                 is close to the liquid-alone reference.

PRESSURE_FEATURES = ['ratio_homog', 'phi2_l']

R_AIR = 287.0530        # gas constant of air, as in the acquisition program [J/(kg K)]

# Po = f Re of fully developed laminar flow in a square duct (Darcy; Shah and London, 1978).
# The friction factor of Churchill (1977) is used with the laminar-equivalent Reynolds
# number of Jones (1976), Re* = (64/Po) Re, continuous over all regimes.
SQUARE_DUCT_PO = 56.91


def two_phase_state(df_experimental, d=8e-3, gas_density='meter'):
    """
    Mass fluxes, quality, homogeneous properties and Reynolds numbers per run.

    Properties at the test-section pressure and temperature (Ptp, Ttp, CoolProp).
    The gas superficial velocity jg is based on the volumetric flow at the orifice
    meter, so the gas mass flux uses the density at the meter (Pg, Tg; ideal gas,
    as in the acquisition program): Gg = Pg / (R_air Tg) jg. With
    ``gas_density='test_section'`` the density at (Ptp, Ttp) is used instead
    (convention up to objective_criteria_8).
    """

    out = properties_df(df_experimental, d=d)

    if gas_density == 'meter':
        rho_meter = df_experimental['Pg'].to_numpy() / (R_AIR * (df_experimental['Tg'].to_numpy() + 273.15))
        out['rhog_meter'] = rho_meter
        out['Gg'] = rho_meter * out['jg'].to_numpy()
        out['g'] = out['Gl'] + out['Gg']
        out['x'] = out['Gg'] / out['g']
    elif gas_density != 'test_section':
        raise ValueError(gas_density)

    out['mu_h'] = 1 / ((1 - out['x']) / out['mul'] + out['x'] / out['mug'])
    out['rho_h'] = 1 / ((1 - out['x']) / out['rhol'] + out['x'] / out['rhog'])
    out['Re_tp'] = out['g'] * d / out['mu_h']
    out['Re_l'] = out['Gl'] * d / out['mul']
    out['Re_g'] = out['Gg'] * d / out['mug']

    return out


def pressure_features(df_experimental, d=8e-3, friction='churchill', gas_density='meter'):
    """
    Pressure-drop features per run from the experimental table (needs run_id,
    Ttp, Ptp, Pg, Tg, jl, jg and dpdz_tp).

    ``friction='churchill'``: Churchill (1977) with the laminar-equivalent Reynolds
    number of a square duct (Jones, 1976). ``friction='piecewise'``: 56.91/Re below
    Re = 2300 and Fang et al. (2011) above (objective_criteria_8, with
    ``gas_density='test_section'``).
    """

    s = two_phase_state(df_experimental, d=d, gas_density=gas_density)

    if friction == 'churchill':
        f_h = friction_churchill(s['Re_tp'].to_numpy(), po=SQUARE_DUCT_PO)
        f_l = friction_churchill(s['Re_l'].to_numpy(), po=SQUARE_DUCT_PO)
    elif friction == 'piecewise':
        f_h = np.array([friction_f(re, laminar_fre=SQUARE_DUCT_PO) for re in s['Re_tp']])
        f_l = np.array([friction_f(re, laminar_fre=SQUARE_DUCT_PO) for re in s['Re_l']])
    else:
        raise ValueError(friction)

    out = pd.DataFrame({
        'run_id': s['run_id'].to_numpy(),
        'dpdz_homog': f_h * s['g'].to_numpy() ** 2 / (2 * d * s['rho_h'].to_numpy()),
        'dpdz_l': f_l * s['Gl'].to_numpy() ** 2 / (2 * d * s['rhol'].to_numpy()),
    }).merge(df_experimental[['run_id', 'dpdz_tp']], on='run_id')

    out['ratio_homog'] = out['dpdz_tp'] / out['dpdz_homog']
    out['phi2_l'] = out['dpdz_tp'] / out['dpdz_l']

    return out[['run_id', 'dpdz_homog', 'dpdz_l'] + PRESSURE_FEATURES]
