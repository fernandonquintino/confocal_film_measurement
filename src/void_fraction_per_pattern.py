import pandas as pd
import numpy as np

# =========================================================
# BUBBLY
# =========================================================
def calculate_void_fraction_bubbly(df, L=8, thick=0.98):
    df = df.copy()
    df['t2'] = df['t2'].fillna(L) # Missing signal interpreted as liquid
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# DISPERSED
# =========================================================
def calculate_void_fraction_dispersed(df, L=8, thick=0.98):
    df = df.copy()
    delta = (df['t1'] - thick).clip(lower=0)
    df['t2'] = df['t2'] + delta
    df['t2'] = df['t2'].fillna(L) # Missing signal interpreted as liquid
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# PLUG
# =========================================================
def calculate_void_fraction_plug(df, L=8, thick=0.98):
    df = df.copy()
    df['t2'] = df['t2'].fillna(L) # Missing signal interpreted as liquid slug
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# SLUG
# =========================================================
def calculate_void_fraction_slug(df, L=8, thick=0.98):
    df = df.copy()
    df['t2'] = df['t2'].fillna(L) # Missing signal interpreted as liquid slug
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# CHURN
# =========================================================
def calculate_void_fraction_churn(df, L=8, thick=0.98):
    df = df.copy()
    delta = (df['t1'] - thick).clip(lower=0)
    df['t2'] = df['t2'] + delta
    vf = 1 - df['t2'] / L # Keep NaNs
    return np.nanmean(vf)

# =========================================================
# STRATIFIED
# =========================================================
def calculate_void_fraction_stratified(df, L=8, thick=0.98):
    df = df.copy()
    delta = (df['t1'] - thick).clip(lower=0)
    df['t2'] = df['t2'] + delta
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# WAVY
# =========================================================
def calculate_void_fraction_wavy(df, L=8, thick=0.98):
    df = df.copy()
    delta = (df['t1'] - thick).clip(lower=0)
    df['t2'] = df['t2'] + delta
    vf = 1 - df['t2'] / L

    return np.nanmean(vf)

# =========================================================
# ANNULAR
# =========================================================

def calculate_void_fraction_annular(df, L=8, thick=0.98):
    df = df.copy()
    delta = (df['t1'] - thick).clip(lower=0)
    df['t2'] = df['t2'] + delta
    vf = 1 - df['t2'] / L # Keep NaNs to avoid artificial film reconstruction

    return np.nanmean(vf)

# =========================================================
# DISPATCHER
# =========================================================

VOID_FRACTION_FUNCTIONS = {
    'churn': calculate_void_fraction_churn,
    'annular': calculate_void_fraction_annular,
    'dispersed': calculate_void_fraction_dispersed,
    'wavy': calculate_void_fraction_wavy,
    'slug': calculate_void_fraction_slug,
    'plug': calculate_void_fraction_plug,
    'stratified': calculate_void_fraction_stratified,
    'bubbly': calculate_void_fraction_bubbly,
}

# =========================================================
# MAIN DISPATCH
# =========================================================

def void_fraction_calculation_function(df, reduced_pattern):

    pattern = str(reduced_pattern).lower()

    if pattern not in VOID_FRACTION_FUNCTIONS:
        raise ValueError(
            f'Unknown pattern: {pattern}'
        )

    return VOID_FRACTION_FUNCTIONS[pattern](df)