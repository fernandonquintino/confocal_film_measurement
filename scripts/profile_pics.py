import os
import sys
from pathlib import Path
import glob

project_root = Path.cwd().parent
sys.path.append(str(project_root))

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns

# # Experimental conditions

silver_path = Path.cwd().parent / 'data/silver_layer/'

silver_confocal_path = Path.cwd().parent / 'data/silver_layer/confocal_runs'
gold_profile_path = Path.cwd().parent / 'data/gold_layer/profile_pics'

gold_profile_path.mkdir(parents=True, exist_ok=True)

df_confocal = pd.read_pickle(silver_path / 'confocal_results.pkl')
df_experimental = pd.read_pickle(silver_path / 'experimental_results.pkl')

# # Selected conditions

missing_threshold = 50

list_runs = df_confocal.loc[df_confocal['t2_missing_percentage']<missing_threshold, 'run_id'].unique()

df_selected = df_experimental.loc[df_experimental['run_id'].isin(list_runs)]

# =========================================================
# Plot settings
# =========================================================

PLOT_STYLE = {
    "ytick.color": "black",
    "xtick.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
}

cm = 1 / 2.54

# Elsevier 5p:
# single column ≈ 8.4 cm
# double column ≈ 17.4 cm

FIG_WIDTH = 10.9 * cm
FIG_HEIGHT = 7.62 * cm

# Optional scaling
FACTOR = 2 / 3

plt.rcParams.update(PLOT_STYLE)

PLOT_CONFIG = {
    # Elsevier 5p two-column
    # single column ≈ 8.4 cm
    "fig_width_cm": 8.4,
    "fig_height_cm": 3.6,

    # Limits
    "ymin": -1.5,
    "ymax": 8.5,

    # Sensor / geometry corrections
    "glass_thickness": 0.96,
    "channel_height": 8.02,

    # t3 correction
    "t3_threshold": 1.3,
    "t3_scale_num": 1.51,
    "t3_scale_den": 1.33,

    # Plot appearance
    "interface_alpha": 0.5,
    "interface_linewidth": 0.8,
    "wall_linewidth": 1.0,

    # Optional time crop
    "max_time": None,   # e.g. 0.3
}

# =========================================================
# Utilities
# =========================================================

def cm_to_inch(value_cm):
    return value_cm / 2.54


def load_confocal_run(run_id, base_path=silver_confocal_path):
    """
    Load confocal run CSV.

    Expected columns:
        t, t1, t2, t3
    """

    file_path = base_path / f'confocal_run_{run_id}.csv'

    df = pd.read_csv(file_path)

    return df

def preprocess_profile(df, config=PLOT_CONFIG):
    """
    Apply legacy corrections and preprocessing.
    """

    df = df.copy()

    # -----------------------------------------------------
    # Remove invalid measurements
    # -----------------------------------------------------
    df = df.loc[~df['t1'].isna()].copy()

    # -----------------------------------------------------
    # Normalize time
    # -----------------------------------------------------
    df['t'] = df['t'] - df['t'].min()

    # -----------------------------------------------------
    # Fill missing t3
    # -----------------------------------------------------
    df.loc[df['t3'].isna(), 't3'] = config['glass_thickness']

    # -----------------------------------------------------
    # Legacy t3 correction
    # -----------------------------------------------------
    mask_t3 = df['t3'] > config['t3_threshold']

    df.loc[mask_t3, 't3'] = (
        df.loc[mask_t3, 't3']
        * config['t3_scale_num']
        / config['t3_scale_den']
    )

    # -----------------------------------------------------
    # Overflow correction
    # -----------------------------------------------------
    mask_overflow = df['t2'] > config['channel_height']

    df.loc[mask_overflow, 't3'] = (
        (df.loc[mask_overflow, 't2'] - config['channel_height'])
        * config['t3_scale_num']
        / config['t3_scale_den']
    )

    df.loc[mask_overflow, 't2'] = config['channel_height']

    # -----------------------------------------------------
    # Bottom wall
    # -----------------------------------------------------
    df['t0'] = -config['glass_thickness']

    # -----------------------------------------------------
    # Optional time crop
    # -----------------------------------------------------
    if config['max_time'] is not None:
        df = df.loc[df['t'] < config['max_time']]

    return df

# =========================================================
# Plot function
# =========================================================

def create_profile_plot(
    df,
    run_id,
    config=PLOT_CONFIG,
    show_xlabel=True,
):
    """
    Create profile plot for a confocal run.
    """

    width = cm_to_inch(config['fig_width_cm'])
    height = cm_to_inch(config['fig_height_cm'])

    fig, ax = plt.subplots(figsize=(width, height))

    # -----------------------------------------------------
    # Walls
    # -----------------------------------------------------
    ax.plot(
        df['t'],
        df['t0'],
        color='k',
        linewidth=config['wall_linewidth']
    )

    ax.plot(
        df['t'],
        df['t1'] - config['glass_thickness'],
        color='k',
        linewidth=config['wall_linewidth']
    )

    ax.plot(
        df['t'],
        config['channel_height'] + 0.1,
        color='k',
        linewidth=config['wall_linewidth']
    )

    # -----------------------------------------------------
    # Interface
    # -----------------------------------------------------
    ax.plot(
        df['t'],
        df['t2'],
        color='blue',
        alpha=config['interface_alpha'],
        linewidth=config['interface_linewidth']
    )

    # -----------------------------------------------------
    # Labels
    # -----------------------------------------------------
    if show_xlabel:
        ax.set_xlabel('Time [s]')
    else:
        ax.set_xlabel(None)

    ax.set_ylabel('Profile [mm]')

    # -----------------------------------------------------
    # Limits
    # -----------------------------------------------------
    ax.set_ylim(config['ymin'], config['ymax'])

    # -----------------------------------------------------
    # Tick formatting
    # -----------------------------------------------------
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    # -----------------------------------------------------
    # Clean style
    # -----------------------------------------------------
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    return fig

# =========================================================
# Save function
# =========================================================

def save_profile_plot(
    fig,
    run_id,
    output_path=gold_profile_path,
    extension='eps',
):
    """
    Save profile plot for publication.
    """

    save_path = output_path / f'profile_run_{run_id}.{extension}'

    fig.savefig(
        save_path,
        format=extension,
        bbox_inches='tight'
    )

    plt.close(fig)

    return save_path

# =========================================================
# Batch processing
# =========================================================

def process_profile_runs(
    run_list,
    config=PLOT_CONFIG,
    input_path=silver_confocal_path,
    output_path=gold_profile_path,
):
    """
    Process multiple confocal runs.
    """

    for run_id in run_list:

        print(f'Processing run {run_id}...')

        # ---------------------------------------------
        # Load
        # ---------------------------------------------
        df = load_confocal_run(
            run_id=run_id,
            base_path=input_path
        )

        # ---------------------------------------------
        # Preprocess
        # ---------------------------------------------
        df = preprocess_profile(
            df=df,
            config=config
        )

        # ---------------------------------------------
        # Plot
        # ---------------------------------------------
        fig = create_profile_plot(
            df=df,
            run_id=run_id,
            config=config
        )

        # ---------------------------------------------
        # Save
        # ---------------------------------------------
        save_profile_plot(
            fig=fig,
            run_id=run_id,
            output_path=output_path
        )

        print(f'Finished run {run_id}')


# =========================================================
# Example
# =========================================================

list_runs = [1, 2, 3, 4]

process_profile_runs(
    run_list=list_runs
)


