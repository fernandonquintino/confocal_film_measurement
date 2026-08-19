import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def map_garimella(df, FIG_WIDTH_CM=8.4, FIG_HEIGHT_CM=7.2):

    CM = 1 / 2.54
    FIG_SIZE = (FIG_WIDTH_CM * CM, FIG_HEIGHT_CM * CM)

    PLOT_STYLE = {

        # Lines
        "line_color": "k",
        "line_width": 0.7,

        # Scatter
        "marker_size": 30,
        "marker_edgecolor": "k",
        "marker_linewidth": 0.5,

        # Grid
        "grid_alpha": 0.5,
    }


    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.grid(alpha=PLOT_STYLE["grid_alpha"], zorder=0)

    jg = pd.Series([0.1049, 1.6286, 0.7571, 2.6015, 4.2504, 17.8180, 13.3123, 6.5641, 16.2489, 98.3638])
    jl = pd.Series([0.7245, 1.2322, 6.0739, 1.2651, 0.0354, 1.1254, 0.0915, 0.0397, 0.0410, 1.0334])

    connections = [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (6, 7), (6, 8),(5, 9),]

    for i, j in connections:

        ax.plot(
            [jg[i], jg[j]],
            [jl[i], jl[j]],
            color=PLOT_STYLE["line_color"],
            linewidth=PLOT_STYLE["line_width"],
            zorder=1
        )

    # =========================================================
    # Scatter points
    # =========================================================

    unique_categories = df["reduced_pattern"].unique()

    markers = ['o', 's', '^', 'v', '*', 'D', 'X', 'P']

    markers_mapping = {
        category: marker
        for category, marker in zip(unique_categories, markers)
    }

    for category in unique_categories:

        subset = df[df["reduced_pattern"] == category]

        ax.scatter(
            subset["jg"],
            subset["jl"],

            s=PLOT_STYLE["marker_size"],

            edgecolor=PLOT_STYLE["marker_edgecolor"],
            linewidth=PLOT_STYLE["marker_linewidth"],

            marker=markers_mapping[category],

            label=category,

            zorder=2
        )

    # =========================================================
    # Labels
    # =========================================================

    ax.text(1.1e-1, 3, "bubbly")
    ax.text(8, 1.5e-2, "stratified")
    ax.text(1.1e-1, 5e-2, "intermittent")
    ax.text(5, 5, "dispersed")
    ax.text(5, 0.5, "wavy")
    ax.text(30, 0.2, "annular")

    # =========================================================
    # Axes
    # =========================================================

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim([0.07, 150])
    ax.set_ylim([0.01, 10])

    ax.set_xlabel(r"$j_g\;$[m/s]")
    ax.set_ylabel(r"$j_l\;$[m/s]")


    ax.legend(title="Flow pattern", bbox_to_anchor=(1, 1), loc="upper left", frameon=False)

    plt.close(fig)

    return fig, ax

def map_garimella2(
    df,
    missing_percentage=50,
    missing_column="t2_missing_percentage",
    FIG_WIDTH_CM=8.4,
    FIG_HEIGHT_CM=7.2,
):
    """
    Same map as ``map_garimella``, but the points are additionally split by
    data quality:

    * marker shape  -> flow pattern (same mapping as ``map_garimella``)
    * marker filling -> missing data

        - filled marker  : ``missing_column`` <  ``missing_percentage``
        - hollow marker  : ``missing_column`` >= ``missing_percentage``

    Encoding the quality with the filling (instead of a second shape or a
    second colour) keeps the flow-pattern legend identical to the one of
    ``map_garimella`` and stays readable in black and white.
    """

    CM = 1 / 2.54
    FIG_SIZE = (FIG_WIDTH_CM * CM, FIG_HEIGHT_CM * CM)

    PLOT_STYLE = {

        # Lines
        "line_color": "k",
        "line_width": 0.7,

        # Scatter
        "marker_size": 14,
        "marker_edgecolor": "k",
        "marker_linewidth": 0.4,
        "marker_linewidth_hollow": 0.7,

        # Grid
        "grid_alpha": 0.5,
    }


    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.grid(alpha=PLOT_STYLE["grid_alpha"], zorder=0)

    jg = pd.Series([0.1049, 1.6286, 0.7571, 2.6015, 4.2504, 17.8180, 13.3123, 6.5641, 16.2489, 98.3638])
    jl = pd.Series([0.7245, 1.2322, 6.0739, 1.2651, 0.0354, 1.1254, 0.0915, 0.0397, 0.0410, 1.0334])

    connections = [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (6, 7), (6, 8),(5, 9),]

    for i, j in connections:

        ax.plot(
            [jg[i], jg[j]],
            [jl[i], jl[j]],
            color=PLOT_STYLE["line_color"],
            linewidth=PLOT_STYLE["line_width"],
            zorder=1
        )

    # =========================================================
    # Scatter points
    # =========================================================

    unique_categories = df["reduced_pattern"].unique()

    markers = ['o', 's', '^', 'v', '*', 'D', 'X', 'P']
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    markers_mapping = {
        category: marker
        for category, marker in zip(unique_categories, markers)
    }

    colors_mapping = {
        category: color
        for category, color in zip(unique_categories, colors)
    }

    is_good = df[missing_column] < missing_percentage

    for category in unique_categories:

        subset = df[df["reduced_pattern"] == category]

        good = subset[is_good.loc[subset.index]]
        poor = subset[~is_good.loc[subset.index]]

        # below the threshold -> filled marker
        ax.scatter(
            good["jg"],
            good["jl"],

            s=PLOT_STYLE["marker_size"],

            color=colors_mapping[category],
            edgecolor=PLOT_STYLE["marker_edgecolor"],
            linewidth=PLOT_STYLE["marker_linewidth"],

            marker=markers_mapping[category],

            zorder=3
        )

        # at or above the threshold -> hollow marker
        ax.scatter(
            poor["jg"],
            poor["jl"],

            s=PLOT_STYLE["marker_size"],

            facecolor="none",
            edgecolor=colors_mapping[category],
            linewidth=PLOT_STYLE["marker_linewidth_hollow"],

            marker=markers_mapping[category],

            zorder=2
        )

    # =========================================================
    # Labels
    # =========================================================

    ax.text(1.1e-1, 3, "bubbly")
    ax.text(8, 1.5e-2, "stratified")
    ax.text(1.1e-1, 5e-2, "intermittent")
    ax.text(5, 5, "dispersed")
    ax.text(5, 0.5, "wavy")
    ax.text(30, 0.2, "annular")

    # =========================================================
    # Axes
    # =========================================================

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim([0.07, 150])
    ax.set_ylim([0.01, 10])

    ax.set_xlabel(r"$j_g\;$[m/s]")
    ax.set_ylabel(r"$j_l\;$[m/s]")

    # =========================================================
    # Legend
    # =========================================================

    pattern_handles = [
        plt.Line2D(
            [], [],
            linestyle="none",
            marker=markers_mapping[category],
            markerfacecolor=colors_mapping[category],
            markeredgecolor=PLOT_STYLE["marker_edgecolor"],
            markeredgewidth=PLOT_STYLE["marker_linewidth"],
            markersize=PLOT_STYLE["marker_size"] ** 0.5,
            label=category,
        )
        for category in unique_categories
    ]

    ax.legend(
        handles=pattern_handles,
        title="Flow pattern",
        bbox_to_anchor=(1, 1),
        loc="upper left",
        frameon=False,
    )

    plt.close(fig)

    return fig, ax
