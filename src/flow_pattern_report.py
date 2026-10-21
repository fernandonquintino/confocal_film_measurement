import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS

CM = 1 / 2.54

FIG_SIZE = (8.4 * CM, 6.8 * CM)        # single column
FIG_SIZE_WIDE = (17.6 * CM, 8.4 * CM)  # double column

VISUAL_PATTERNS = ['bubbly', 'plug', 'slug', 'dispersed', 'churn', 'wavy', 'stratified', 'annular']


def pattern_styles(map_categories):
    """
    Marker and colour per pattern, in the same order as the flow pattern
    maps (``map_garimella``): the i-th category takes the i-th marker and
    the i-th colour of the default cycle.
    """

    markers = ['o', 's', '^', 'v', '*', 'D', 'X', 'P']
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    return {
        category: {'marker': marker, 'color': color}
        for category, marker, color in zip(map_categories, markers, colors)
    }


def _pattern_legend(ax, styles, patterns, title='Visual pattern', loc=None):
    """
    Legend of the flow patterns. ``loc=None`` places it outside the axes, on
    the right; any matplotlib location (e.g. 'best') places it inside.
    """

    handles = [
        Line2D([], [], linestyle='none', marker=styles[p]['marker'], markerfacecolor=styles[p]['color'],
               markeredgecolor='k', markeredgewidth=0.4, markersize=4, label=p)
        for p in patterns
    ]

    if loc is None:
        ax.legend(handles=handles, title=title, bbox_to_anchor=(1, 1), loc='upper left',
                  frameon=False, fontsize=7, title_fontsize=7)
    else:
        legend = ax.legend(handles=handles, title=title, loc=loc, frameon=True, framealpha=0.9, edgecolor='0.6',
                           fancybox=False, fontsize=7, title_fontsize=7, borderpad=0.3, labelspacing=0.25,
                           handletextpad=0.2)
        legend.get_frame().set_linewidth(0.5)

# =========================================================
# Criterion figures
# =========================================================

def criterion_scatter(df, x, y, xlabel, ylabel, styles, xlines=(), ylines=(), xscale='linear',
                      pattern_column='reduced_pattern', patterns=VISUAL_PATTERNS, segments=(),
                      xlim=None, ylim=None, yscale='linear', label_column=None,
                      figsize=None, legend_loc=None):
    """
    Scatter of two run features coloured by pattern. ``xlines``/``ylines``
    draw full dashed threshold lines; ``segments`` draws dashed segments
    given as ((x0, x1), (y0, y1)); ``label_column`` annotates each point
    (e.g. with the run id). ``figsize`` overrides the default size and
    ``legend_loc`` places the legend inside the axes (see ``_pattern_legend``).
    """

    fig, ax = plt.subplots(figsize=figsize or FIG_SIZE)

    present = [p for p in patterns if (df[pattern_column] == p).any()]

    for pattern in present:

        subset = df[df[pattern_column] == pattern]

        ax.scatter(subset[x], subset[y], s=16, marker=styles[pattern]['marker'],
                   color=styles[pattern]['color'], edgecolor='k', linewidth=0.4, zorder=2)

    for xline in xlines:
        ax.axvline(xline, color='k', linestyle='--', linewidth=0.7, zorder=1)

    for yline in ylines:
        ax.axhline(yline, color='k', linestyle='--', linewidth=0.7, zorder=1)

    for xs, ys in segments:
        ax.plot(xs, ys, color='k', linestyle='--', linewidth=0.7, zorder=1)

    if label_column is not None:
        for _, row in df[df[pattern_column].isin(present)].iterrows():
            ax.annotate(str(row[label_column]), (row[x], row[y]), xytext=(3, 2),
                        textcoords='offset points', fontsize=5, zorder=3)

    if xscale == 'symlog':
        ax.set_xscale('symlog', linthresh=0.1)
        ax.set_xlim(-0.01, 1.5 * max(df[x].max(), 0.2))

    if yscale != 'linear':
        ax.set_yscale(yscale)

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3, zorder=0)

    _pattern_legend(ax, styles, present, loc=legend_loc)

    fig.tight_layout()
    plt.close(fig)

    return fig


def confusion_heatmap(confusion, xlabel='Objective pattern', ylabel='Visual pattern', figsize=None):

    fig, ax = plt.subplots(figsize=figsize or (8.4 * CM, 8.4 * CM))

    sns.heatmap(confusion, annot=True, fmt='d', cmap='Greys', cbar=False, linewidths=0.5,
                linecolor='white', annot_kws={'fontsize': 7}, ax=ax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='x', rotation=30)
    ax.tick_params(axis='y', rotation=0)

    fig.tight_layout()
    plt.close(fig)

    return fig


def sweep_figure(df_sweep, defaults, labels, log_parameters=('f_bridge_hz',),
                 width=17.6 * CM, row_height=4.2 * CM, label_size=7, tick_size=6,
                 ylabel=r'Balanced accuracy [\%]'):
    """
    Balanced accuracy vs. each parameter (one panel per parameter, default
    value marked with a dashed line). ``width`` and ``row_height`` are in
    inches.
    """

    parameters = list(dict.fromkeys(df_sweep['parameter']))
    ncols = 4
    nrows = int(np.ceil(len(parameters) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(width, row_height * nrows), sharey=True, squeeze=False)

    for ax, key in zip(axes.flat, parameters):

        group = df_sweep[df_sweep['parameter'] == key]

        ax.plot(group['value'], 100 * group['balanced_accuracy'], color='#4C72B0', linewidth=1.2, marker='o', markersize=3)
        ax.axvline(defaults[key], color='k', linestyle='--', linewidth=0.7)

        if key in log_parameters:
            ax.set_xscale('log')

        ax.set_xlabel(labels[key], fontsize=label_size)
        ax.tick_params(labelsize=tick_size)
        ax.grid(alpha=0.3)

    for ax in axes.flat[len(parameters):]:
        ax.set_visible(False)

    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=label_size)

    fig.tight_layout()
    plt.close(fig)

    return fig


def bar_figure(values, xlabel, color='#4C72B0', highlight=None, highlight_color='#C44E52',
               legend_labels=('missing-data feature', 'detected-thickness feature'), groups=None):
    """
    Horizontal bar chart of a Series (index = LaTeX-safe labels). Bars whose
    index is in ``highlight`` get ``highlight_color``; ``legend_labels``
    names the highlighted and the regular bars. ``groups`` (legend label ->
    (members, colour)) adds further coloured groups.
    """

    values = values.sort_values()
    highlight = set(highlight or [])
    groups = groups or {}

    fig, ax = plt.subplots(figsize=(8.4 * CM, max(4.0, 0.42 * len(values)) * CM))

    def bar_color(label):
        for members, group_color in groups.values():
            if label in members:
                return group_color
        return highlight_color if label in highlight else color

    colors = [bar_color(label) for label in values.index]

    ax.barh(range(len(values)), values.to_numpy(), color=colors, height=0.7)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(values.index, fontsize=7)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis='x', labelsize=7)
    ax.grid(axis='x', alpha=0.3)

    if highlight or groups:
        handles = [
            plt.Rectangle((0, 0), 1, 1, color=highlight_color, label=legend_labels[0]),
            plt.Rectangle((0, 0), 1, 1, color=color, label=legend_labels[1]),
        ] if highlight else []

        handles += [
            plt.Rectangle((0, 0), 1, 1, color=group_color, label=label)
            for label, (_, group_color) in groups.items()
        ]

        ax.legend(handles=handles, loc='lower right', fontsize=6, frameon=False)

    fig.tight_layout()
    plt.close(fig)

    return fig

# =========================================================
# Clustering figures
# =========================================================

def mds_figure(distance, patterns, styles, title_patterns=VISUAL_PATTERNS, seed=0):
    """
    2D metric MDS embedding of a condensed distance matrix, coloured by
    visual pattern.
    """

    embedding = MDS(n_components=2, metric='precomputed', n_init=4, init='random',
                    random_state=seed).fit_transform(squareform(distance))

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    patterns = np.asarray(patterns)
    present = [p for p in title_patterns if (patterns == p).any()]

    for pattern in present:

        sel = patterns == pattern

        ax.scatter(embedding[sel, 0], embedding[sel, 1], s=16, marker=styles[pattern]['marker'],
                   color=styles[pattern]['color'], edgecolor='k', linewidth=0.4, zorder=2)

    ax.set_xlabel('MDS 1 [-]')
    ax.set_ylabel('MDS 2 [-]')
    ax.grid(alpha=0.3, zorder=0)

    _pattern_legend(ax, styles, present)

    fig.tight_layout()
    plt.close(fig)

    return fig


def dendrogram_figure(distance, patterns, styles, k, method='average', title_patterns=VISUAL_PATTERNS):
    """
    Dendrogram of a condensed distance matrix with a colour strip showing
    the visual pattern of each leaf, cut into ``k`` clusters.
    """

    Z = linkage(distance, method=method)
    cut = Z[-(k - 1), 2] if k > 1 else Z[-1, 2] * 1.01

    fig, (ax, strip) = plt.subplots(
        2, 1, figsize=FIG_SIZE_WIDE, gridspec_kw={'height_ratios': [6, 0.6], 'hspace': 0.02}, sharex=False
    )

    tree = dendrogram(Z, ax=ax, no_labels=True, color_threshold=cut, above_threshold_color='0.5')

    ax.axhline(cut, color='k', linestyle='--', linewidth=0.7)
    ax.set_ylabel('Distance [-]')
    ax.set_xticks([])

    patterns = np.asarray(patterns)[tree['leaves']]
    colors = [styles[p]['color'] if p in styles else '0.8' for p in patterns]

    for i, color in enumerate(colors):
        strip.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, linewidth=0))

    strip.set_xlim(0, len(colors))
    strip.set_ylim(0, 1)
    strip.axis('off')

    present = [p for p in title_patterns if (patterns == p).any()]
    _pattern_legend(ax, styles, present)

    plt.close(fig)

    return fig

# =========================================================
# Profile grid
# =========================================================

def _densest_window(df, window):
    """
    Start time of the ``window``-long interval with the most detected t2.
    """

    time = df['time'].to_numpy()
    detected = df['t2'].notna().to_numpy().astype(int)

    if len(time) == 0 or time[-1] <= window:
        return 0.0

    cumulative = np.concatenate([[0], np.cumsum(detected)])
    ends = np.searchsorted(time, time + window)
    counts = cumulative[ends] - cumulative[np.arange(len(time))]

    valid = time + window <= time[-1]

    return float(time[np.argmax(np.where(valid, counts, -1))])


def profile_grid_figure(panels, film_estimate_fn, window=1.0, ncols=3, y_max=3.0):
    """
    Grid of film profiles. ``panels`` is a list of (df, title). Each panel
    shows, over the ``window``-long interval with the most detected samples,
    the detected t2 (blue dots), the film below the detection limit
    estimated from t1 (orange dots, from ``film_estimate_fn``) and the wall.
    Readings above ``y_max`` are clipped.
    """

    nrows = int(np.ceil(len(panels) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(17.6 * CM, 3.8 * nrows * CM), sharey=True, squeeze=False)

    for ax, (df, title) in zip(axes.flat, panels):

        if 'mixer_wave' in df.columns:
            df = df.loc[~df['mixer_wave'].astype(bool)].reset_index(drop=True)

        start = _densest_window(df, window)
        df = df[(df['time'] >= start) & (df['time'] <= start + window)].copy()
        df['time'] = df['time'] - start

        estimate = film_estimate_fn(df)
        merged = df['t2'].isna().to_numpy() & ~np.isnan(estimate)

        ax.plot(df['time'], np.zeros(len(df)), color='k', linewidth=0.6)
        ax.plot(df['time'], df['t2'], '.', color='#4C72B0', markersize=1.2)
        ax.plot(df['time'].to_numpy()[merged], estimate[merged], '.', color='#DD8452', markersize=1.2)

        ax.set_title(title, fontsize=5.5, loc='left')
        ax.tick_params(labelsize=5)
        ax.set_ylim(-0.1, y_max)

    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)

    for ax in axes[:, 0]:
        ax.set_ylabel('Profile [mm]', fontsize=6)

    for ax in axes[-1, :]:
        ax.set_xlabel('Time [s]', fontsize=6)

    handles = [
        Line2D([], [], linestyle='none', marker='.', color='#4C72B0', label='detected $t_2$'),
        Line2D([], [], linestyle='none', marker='.', color='#DD8452', label='film below detection (from $t_1$)'),
    ]
    fig.legend(handles=handles, loc='upper center', ncol=2, fontsize=6, frameon=False, bbox_to_anchor=(0.5, 1.02))

    fig.tight_layout()
    plt.close(fig)

    return fig

# =========================================================
# Mixer wave figure
# =========================================================

def profile_comparison_figure(df_before, df_after, t_max=2.0, glass_thickness=0.96):
    """
    t2 profile before (processed) and after (processed_v2) the mixer wave
    removal, with the removed samples highlighted.
    """

    fig, axes = plt.subplots(2, 1, figsize=(17.6 * CM, 7.0 * CM), sharex=True, sharey=True)

    before = df_before[df_before['time'] <= t_max]
    after = df_after[df_after['time'] <= t_max]
    removed = after['mixer_wave'].astype(bool).to_numpy()

    for ax, df, title in [(axes[0], before, 'processed'), (axes[1], after, 'processed\\_v2')]:

        ax.plot(df['time'], df['t1'] - glass_thickness, color='k', linewidth=0.6)
        ax.plot(df['time'], df['t2'], '.', color='#4C72B0', markersize=1.5)
        ax.set_ylabel('Profile [mm]', fontsize=7)
        ax.set_title(title, fontsize=7, loc='left')
        ax.tick_params(labelsize=6)

    axes[0].plot(before['time'].to_numpy()[removed], before['t2'].to_numpy()[removed], '.',
                 color='#C44E52', markersize=2.5, label='removed (mixer wave)')
    axes[0].legend(fontsize=6, frameon=False, loc='upper right')

    axes[1].set_xlabel('Time [s]', fontsize=7)
    axes[1].set_ylim(-0.2, 8.4)

    fig.tight_layout()
    plt.close(fig)

    return fig
