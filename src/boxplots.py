import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

def boxplot_plot(df, column, label, unit, FIG_WIDTH_CM = 8.4, FIG_HEIGHT_CM = 6.8, patterns=['bubbly', 'slug', 'plug', 'wavy', 'dispersed']):

    # Order of patterns

    CM = 1 / 2.54
    FIG_SIZE = (FIG_WIDTH_CM * CM, FIG_HEIGHT_CM * CM)

    # Extract grouped data
    grouped_data = [
        df.loc[
            df['reduced_pattern'] == pattern,
            column
        ].dropna()
        for pattern in patterns
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    # More serious academic colors
    colors = [
        '#4C72B0',  # blue
        '#55A868',  # green
        '#C44E52',  # red
        '#8172B2',  # purple
        '#CCB974',   # yellow-brown
        '#937860',   # brown
        '#8C8C8C',  # gray
        '#DD8452',  # orange

    ]

    # Boxplot
    bp = ax.boxplot(
        grouped_data,
        patch_artist=True,
        widths=0.6,
        showfliers=False,
    )

    # Fill colors
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
        patch.set_linewidth(1.1)

    # Style other elements
    for median in bp['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)

    for whisker in bp['whiskers']:
        whisker.set_linewidth(1.0)

    for cap in bp['caps']:
        cap.set_linewidth(1.0)

    # Labels
    ax.set_xticks(range(1, len(patterns) + 1))
    ax.set_xticklabels(patterns, rotation=20, ha='right')

    ax.set_xlabel(r'Flow pattern')
    ax.set_ylabel(f'{label} [{unit}]')

    # Cleaner style
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # Subtle horizontal grid
    ax.grid(
        axis='y',
        linestyle='--',
        linewidth=0.6,
        alpha=0.35
    )

    # Layout adjustment
    fig.tight_layout()

    plt.close()

    return fig, ax