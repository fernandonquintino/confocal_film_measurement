import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull

# =========================================================
# Coleman & Garimella (1999) transition lines
# =========================================================

GARIMELLA_JG = pd.Series([0.1049, 1.6286, 0.7571, 2.6015, 4.2504, 17.8180, 13.3123, 6.5641, 16.2489, 98.3638])
GARIMELLA_JL = pd.Series([0.7245, 1.2322, 6.0739, 1.2651, 0.0354, 1.1254, 0.0915, 0.0397, 0.0410, 1.0334])
GARIMELLA_CONNECTIONS = [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (6, 7), (6, 8), (5, 9)]


def _draw_garimella_lines(ax, style):

    for i, j in GARIMELLA_CONNECTIONS:

        ax.plot(
            [GARIMELLA_JG[i], GARIMELLA_JG[j]],
            [GARIMELLA_JL[i], GARIMELLA_JL[j]],
            color=style["line_color"],
            linewidth=style["line_width"],
            zorder=1
        )


def _draw_garimella_labels(ax):

    ax.text(1.1e-1, 3, "bubbly")
    ax.text(8, 1.5e-2, "stratified")
    ax.text(1.1e-1, 5e-2, "intermittent")
    ax.text(5, 5, "dispersed")
    ax.text(5, 0.5, "wavy")
    ax.text(30, 0.2, "annular")


def _set_map_axes(ax):

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim([0.07, 150])
    ax.set_ylim([0.01, 10])

    ax.set_xlabel(r"$j_g\;$[m/s]")
    ax.set_ylabel(r"$j_l\;$[m/s]")


def map_garimella(df, FIG_WIDTH_CM=8.4, FIG_HEIGHT_CM=7.2, style_order=None):
    """
    ``style_order`` (optional) fixes which marker/colour each category gets:
    the i-th category of the list takes the i-th marker and the i-th colour
    of the cycle. Without it, styles follow the order of appearance in ``df``.
    """

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

    _draw_garimella_lines(ax, PLOT_STYLE)

    # =========================================================
    # Scatter points
    # =========================================================

    unique_categories = df["reduced_pattern"].unique()

    markers = ['o', 's', '^', 'v', '*', 'D', 'X', 'P']

    style_categories = unique_categories if style_order is None else style_order

    markers_mapping = {
        category: marker
        for category, marker in zip(style_categories, markers)
    }

    colors_mapping = {
        category: color
        for category, color in zip(
            style_categories,
            plt.rcParams["axes.prop_cycle"].by_key()["color"],
        )
    }

    for category in unique_categories:

        subset = df[df["reduced_pattern"] == category]

        color_kwargs = {} if style_order is None else {"color": colors_mapping[category]}

        ax.scatter(
            subset["jg"],
            subset["jl"],

            s=PLOT_STYLE["marker_size"],

            edgecolor=PLOT_STYLE["marker_edgecolor"],
            linewidth=PLOT_STYLE["marker_linewidth"],

            marker=markers_mapping[category],

            label=category,

            zorder=2,

            **color_kwargs
        )

    # =========================================================
    # Labels
    # =========================================================

    _draw_garimella_labels(ax)

    # =========================================================
    # Axes
    # =========================================================

    _set_map_axes(ax)


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

    _draw_garimella_lines(ax, PLOT_STYLE)

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

    _draw_garimella_labels(ax)

    # =========================================================
    # Axes
    # =========================================================

    _set_map_axes(ax)

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

def map_garimella_region(
    df,
    pad=0.06,
    label="evaluated\nregion",
    label_xy=(0.4, 0.3),
    FIG_WIDTH_CM=8.4,
    FIG_HEIGHT_CM=7.2,
):
    """
    Same map as ``map_garimella``, but instead of plotting each run the whole
    tested envelope is shown as a single shaded region.

    The region is the convex hull of the (log10 jg, log10 jl) points of
    ``df``. Every point is first replaced by a small circle of radius ``pad``
    (in decades) so the runs on the edge fall inside the hull and not on its
    outline.
    """

    CM = 1 / 2.54
    FIG_SIZE = (FIG_WIDTH_CM * CM, FIG_HEIGHT_CM * CM)

    PLOT_STYLE = {

        # Lines
        "line_color": "k",
        "line_width": 0.7,

        # Region
        "region_facecolor": "0.6",
        "region_alpha": 0.3,
        "region_edgecolor": "0.25",
        "region_linewidth": 0.8,

        # Grid
        "grid_alpha": 0.5,
    }

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.grid(alpha=PLOT_STYLE["grid_alpha"], zorder=0)

    _draw_garimella_lines(ax, PLOT_STYLE)

    # =========================================================
    # Evaluated region (padded convex hull in log space)
    # =========================================================

    log_points = np.log10(df[["jg", "jl"]].dropna().to_numpy())

    angles = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    circle = pad * np.column_stack([np.cos(angles), np.sin(angles)])

    padded = (log_points[:, None, :] + circle[None, :, :]).reshape(-1, 2)

    hull = ConvexHull(padded)
    vertices = 10 ** padded[hull.vertices]

    ax.add_patch(
        Polygon(
            vertices,
            closed=True,
            facecolor=PLOT_STYLE["region_facecolor"],
            alpha=PLOT_STYLE["region_alpha"],
            edgecolor="none",
            zorder=0.5,
        )
    )

    ax.add_patch(
        Polygon(
            vertices,
            closed=True,
            facecolor="none",
            edgecolor=PLOT_STYLE["region_edgecolor"],
            linewidth=PLOT_STYLE["region_linewidth"],
            linestyle="--",
            zorder=1.5,
        )
    )

    ax.text(
        label_xy[0],
        label_xy[1],
        label,
        ha="center",
        va="center",
        multialignment="center",
        fontsize="small",
        zorder=3,
    )

    # =========================================================
    # Labels
    # =========================================================

    _draw_garimella_labels(ax)

    # =========================================================
    # Axes
    # =========================================================

    _set_map_axes(ax)

    plt.close(fig)

    return fig, ax
