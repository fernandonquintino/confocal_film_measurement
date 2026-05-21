import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_percentage_error, r2_score

def eta_n(y_true, y_pred, n):
    eta = y_true[((y_true - y_pred).abs()/y_true).abs() <= n].shape[0]/y_true.shape[0]
    return eta

def find_closest_category(value, categories):
    return min(categories, key=lambda x: abs(x - value))

def plot_category_data(ax, df, categories, j_cat, j_x, positions, condition, remove_index, xlabel, ylabel, legend_title):
    """
    Plot category data on a given axis.

    Parameters:
        ax: matplotlib axis
        df: DataFrame containing the data
        categories: list of categories to plot
        positions: list of category positions
        condition: lambda function for filtering data
        remove_index: indices to exclude
        xlabel, ylabel: labels for the axes
        legend_title: title for the legend
    """
    palette = sns.color_palette(n_colors=len(categories))
    # marker = ['o', 's', '^', 'v', '*', 'D', 'X', 'P']

    for i, category in enumerate(categories[positions]):
    # for category in categories[positions]:
        category_data = df[
            (df[j_cat] == category) & #(df['jg_cat'] == category) &
            condition(df) &
            (~df.index.isin(remove_index))
        ]
        sns.lineplot(
            x=category_data[j_x], #x=category_data['jl'],
            y=category_data['dpdz_exp_round'],
            marker='o',# marker[i],#'o',
            label=f'{category}',
            ax=ax
        )
        # sns.lineplot(
        #     x=category_data[j_x], #x=category_data['jl'],
        #     y=category_data['muller_heck']/1000,
        #     ax=ax,
        #     color=palette[i],
        #     linestyle='--',
        # )
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.legend(title=legend_title)
    ax.grid(True)

def table_metrics(dfz, list_authors):
    df_metrics = pd.DataFrame({'MAPE': [], 'eta_10': [], 'eta_30': [], 'eta_50': []})
    y_true = dfz['dpdz_exp']
    for author in list_authors:
        y_pred = dfz[author]
        mape = mean_absolute_percentage_error(y_true, y_pred)
        eta10 = eta_n(y_true, y_pred, 0.1)
        eta30 = eta_n(y_true, y_pred, 0.3)
        eta50 = eta_n(y_true, y_pred, 0.5)
        y_mean = y_true.mean()
        df_author = pd.DataFrame({'MAPE': [mape], 'eta_10': [eta10], 'eta_30': [eta30], 'eta_50': [eta50]}, index=[author])
        df_metrics = pd.concat([df_metrics, df_author])
        # df_metrics.sort_values(by='MAPE', inplace=True)
    return(df_metrics)

def table_metrics2(dfz, list_authors):
    df_metrics = pd.DataFrame({'MAPE': [], 'eta_10': [], 'eta_30': [], 'eta_50': []})
    y_true = dfz['dpdz_exp']
    for author in list_authors:
        y_pred = dfz[author]
        mape = mean_absolute_percentage_error(y_true, y_pred)
        eta10 = eta_n(y_true, y_pred, 0.1)
        eta30 = eta_n(y_true, y_pred, 0.3)
        eta50 = eta_n(y_true, y_pred, 0.5)
        y_mean = y_true.mean()
        df_author = pd.DataFrame({'MAPE': [mape], 'eta_10': [eta10], 'eta_30': [eta30], 'eta_50': [eta50]}, index=[author])
        df_metrics = pd.concat([df_metrics, df_author])
    df_metrics = df_metrics.loc[['mc_addams', 'dukler', 'cicchitti', 'chisholm', 'friedel', 'beatti_walley', 'muller_heck', 'lin', 'michima_hibiki', 'chen', 'sun_mishima', 'tibirica', 'nie']]
    df_metrics = df_metrics[['MAPE', 'eta_30']]
    return(df_metrics)

def re_corr_plot(df, author, Re=2000, cm=1/2.56, factor=2/3, bin_edges=[0, 12, 30, 45, 60, 80, np.inf], bin_labels=["0-12", "12-30", "30-45", "45-60", "60-80", "80+"]):

    df['mpe'] = np.abs((df['dpdz_exp']-df[author])/df['dpdz_exp'])*100
    df['mpe_cat'] = pd.cut(df['mpe'], bins=bin_edges, labels=bin_labels, include_lowest=True)

    unique_categories = df['mpe_cat'].sort_values().unique()
    markers_mapping = {category: marker for category, marker in zip(unique_categories, ['o', 's', '^', 'v', '*', 'D', 'X', 'P'])}

    fig, ax = plt.subplots(figsize=(10.9*cm*factor, 7.62*cm))
    ax.plot([Re, Re], [1*10**2, 8*10**4], linestyle='-', alpha=1, color='r', linewidth=0.8, zorder=1)
    ax.plot([1*10**2, 8*10**4], [Re, Re], linestyle='-', alpha=1, color='r', linewidth=0.8, zorder=1)

    plt.text(1.5*10**2, 2*10**2, r'll', size=12)
    plt.text(1.5*10**2, 5*10**4, r'lt', size=12)
    plt.text(5*10**4, 2*10**2, r'tl', size=12)
    plt.text(5*10**4, 5*10**4, r'tt', size=12)

    for category in unique_categories:
        subset = df.sort_values(by='mpe_cat')[df['mpe_cat'] == category]
        ax.scatter(x=subset['Rel'], y=subset['Reg'], s=20, edgecolor='k', linewidth=0.8, marker=markers_mapping[category], label=category, zorder=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlim([1*10**2, 8*10**4])
    plt.ylim([1*10**2, 8*10**4])

    ax.legend(title=r'MAPE$\;$[\%]', bbox_to_anchor=(1, 1.), frameon=False, borderaxespad=0)
    ax.set(xlabel=r'$Re_l$', ylabel='')
    ax.tick_params(axis="both", which="both", direction="in")

    # plt.tight_layout()
    plt.close(fig)  # Prevent redundant display
    return fig

def re_corr_plot_no_legend(df, author, Re=2000, cm=1/2.56, factor=2/3, bin_edges=[0, 12, 30, 45, 60, 80, np.inf], bin_labels=["0-12", "12-30", "30-45", "45-60", "60-80", "80+"]):

    df['mpe'] = np.abs((df['dpdz_exp']-df[author])/df['dpdz_exp'])*100
    df['mpe_cat'] = pd.cut(df['mpe'], bins=bin_edges, labels=bin_labels, include_lowest=True)

    unique_categories = df['mpe_cat'].sort_values().unique()
    markers_mapping = {category: marker for category, marker in zip(unique_categories, ['o', 's', '^', 'v', '*', 'D', 'X', 'P'])}

    fig, ax = plt.subplots(figsize=(10.9*cm*factor, 7.62*cm))
    ax.plot([Re, Re], [1*10**2, 8*10**4], linestyle='-', alpha=1, color='r', linewidth=0.8, zorder=1)
    ax.plot([1*10**2, 8*10**4], [Re, Re], linestyle='-', alpha=1, color='r', linewidth=0.8, zorder=1)

    plt.text(1.5*10**2, 2*10**2, r'll', size=12)
    plt.text(1.5*10**2, 5*10**4, r'lt', size=12)
    plt.text(5*10**4, 2*10**2, r'tl', size=12)
    plt.text(5*10**4, 5*10**4, r'tt', size=12)

    for category in unique_categories:
        subset = df.sort_values(by='mpe_cat')[df['mpe_cat'] == category]
        ax.scatter(x=subset['Rel'], y=subset['Reg'], s=20, edgecolor='k', linewidth=0.8, marker=markers_mapping[category], label=category, zorder=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlim([1*10**2, 8*10**4])
    plt.ylim([1*10**2, 8*10**4])

    # ax.legend(title=r'MAPE$\;$[\%]', bbox_to_anchor=(1., 1.))
    ax.set(xlabel=r'$Re_l$', ylabel=r'$Re_g$')
    ax.tick_params(axis="both", which="both", direction="in")

    plt.tight_layout()
    plt.close(fig)  # Prevent redundant display
    return fig

def tube_profile(et, cm=1/2.56, factor=2/3, time_limit=6, ylim=(-0.5, 10.5), color='k', shade_color='LightBlue'):
    """
    Plots the tube profile with shaded regions and annotations for mean differences.

    Parameters:
        et (pd.DataFrame): Input DataFrame with columns ['time', 't0', 't1', 't2', 't3'].
        cm (float): Scaling factor for figure dimensions.
        factor (float): Adjustment factor for aspect ratio.
        time_limit (float): Maximum time to consider for plotting.
        ylim (tuple): y-axis limits for the plot.
        color (str): Line color.
        shade_color (str): Shading color between regions.

    Returns:
        fig (matplotlib.figure.Figure): The generated plot figure.
    """
    # Filter and modify data
    df_empty = et.loc[(et['t2']>7) & (et['t3'] > 0.6)].copy()
    df_empty['t0'] = 0
    df_empty['t2'] = df_empty['t2'] + df_empty['t1']
    df_empty['t3'] = df_empty['t3'] + df_empty['t2']
    df_empty = df_empty.loc[df_empty['time'] < time_limit]

    # Calculate means
    mean_t1_t0 = (df_empty['t1'] - df_empty['t0']).mean()
    mean_t2_t1 = (df_empty['t2'] - df_empty['t1']).mean()
    mean_t3_t2 = (df_empty['t3'] - df_empty['t2']).mean()

    # Create the plot
    fig = plt.figure(figsize=(10.9 * cm * factor, 7.62 * cm))

    # Plot lines
    plt.plot(df_empty['time'], df_empty['t0'], label='t0', color=color)
    plt.plot(df_empty['time'], df_empty['t1'], label='t1', color=color)
    plt.plot(df_empty['time'], df_empty['t2'], label='t2', color=color)
    plt.plot(df_empty['time'], df_empty['t3'], label='t3', color=color)

    # Shade areas between t0 and t1, t2 and t3
    plt.fill_between(df_empty['time'], df_empty['t0'], df_empty['t1'], color=shade_color, alpha=0.5)
    plt.fill_between(df_empty['time'], df_empty['t2'], df_empty['t3'], color=shade_color, alpha=0.5)

    # Customize axis
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Tube profile [mm]')
    ax.set_ylim(ylim)

    # Add annotations for means
    plt.annotate(f"Mean: {mean_t1_t0:.3f}", xy=(0.05, 1.3), fontsize=10)
    plt.annotate(f"Mean: {mean_t2_t1:.3f}", xy=(0.05, 5), fontsize=10)
    plt.annotate(f"Mean: {mean_t3_t2:.3f}", xy=(0.05, 8), fontsize=10)

    # Return the figure
    plt.tight_layout()
    plt.close(fig)  # Prevent redundant display
    return fig

def list_selected_dataframes(z, df):
    z_select = []

    for _, row in df[['df2', 'number2']].iterrows():
        outer_index = row['df2']  # Index for the outer list
        inner_index = row['number2']  # Index for the inner list

        # Check if the indices are valid
        if 0 <= outer_index < len(z) and 0 <= inner_index < len(z[outer_index]):
            # Append the corresponding DataFrame to z_plug
            z_select.append(z[outer_index][inner_index])
    return z_select

def profile_plot(z, i, df, time=1, cm=1/2.56, factor_h=1.5, factor_v=0.7):
    df3 = z.copy()
    df3 = df3.loc[~df3['t1'].isna()]
    df3['time'] = df3['time'] - df3['time'].min()
    df3.loc[df3['t3'].isna(), 't3'] = 0.96
    df3.loc[df3['t3']>1.3, 't3'] = df3.loc[df3['t3']>1.3, 't3']*1.51/1.33
    df3.loc[df3['t2']>8.02, 't3'] = (df3['t2']-8.02)*1.51/1.33
    df3.loc[df3['t2']>8.02, 't2'] = 8.02
    df3['t0'] = 0-0.96
    df3 = df3.loc[df3['time']<time]

    # Create the plot
    fig = plt.figure(figsize=(10.9*cm*factor_h, 7.62*cm*factor_v))

    # Plot lines
    sns.lineplot(x=df3['time'], y=df3['t0'], color='k')
    sns.lineplot(x=df3['time'], y=df3['t1']-0.96, color='k')
    sns.lineplot(x=df3['time'], y=df3['t2'], color='blue', marker='o', markersize=3, alpha=0.5)#'#89CFF0'
    sns.lineplot(x=df3['time'], y=8.02+0.1, color='k')
    sns.lineplot(x=df3['time'], y=8.98, color='k')

    # Customize axis
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    if i==0:
        ax.set_xlabel('Time [s]')
    else:
        ax.set_xlabel(None)
    ax.set_ylabel('Tube profile [mm]')
    ax.set_ylim(-1.5, 9.5)
    # ax.set_title(f"{i}: $j_l$: {df['jl'][i]}, $j_g$: {df['jg'][i]}")

    # Return the figure
    plt.tight_layout()
    plt.close(fig)  # Prevent redundant display
    return fig

def kde_plot(z, i, df, cm=1/2.56, factor3=1):
    df3 = z.copy()
    df3 = df3.loc[~df3['t1'].isna()]
    df3['time'] = df3['time'] - df3['time'].min()
    df3.loc[df3['t2']>8.02, 't2'] = 8.02

    # Create the plot
    fig = plt.figure(figsize=(10.9 * cm * factor3, 7.62 * cm))

    # Plot lines
    sns.kdeplot(x=df3['t2'], color='k')

    # Customize axis
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.set_xlabel('Film thickness [mm]')
    # ax.set_title(f"{i}: $j_l$: {df['jl'][i]}, $j_g$: {df['jg'][i]}")

    # Return the figure
    plt.tight_layout()
    plt.close(fig)  # Prevent redundant display
    return fig


def kde_combined(z, df, cm=1/2.56, factor3=1):
    fig = plt.figure(figsize=(10.9 * cm * factor3, 7.62 * cm))
    ax = plt.gca()
    palette = sns.color_palette("deep", 7)

    for i in range(len(z)):
        df3 = z[i].copy()
        df3 = df3.loc[~df3['t1'].isna()]
        df3['time'] = df3['time'] - df3['time'].min()
        df3.loc[df3['t2']>8.02, 't2'] = 8.02
        sns.kdeplot(x=df3['t2'], ax=ax, color=palette[i],
                    label=f"$j_l$: {df.iloc[i]['jl']} m/s, $j_g$: {df.iloc[i]['jg']} m/s")


    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.set_xlabel('Film thickness [mm]')
    # ax.legend(bbox_to_anchor=(1.05, 1))
    ax.legend(    loc='upper left',
    bbox_to_anchor=(1.05, 1),  # Outside the plot
    borderaxespad=0, frameon=False)

    # plt.tight_layout()
    plt.close(fig)

    return fig


def kde_combined2(z, df, cm=1/2.56, factor3=1):
    fig = plt.figure(figsize=(10.9 * cm * factor3, 7.62 * cm))
    ax = plt.gca()
    palette = sns.color_palette("deep", 7)

    for i in range(len(z)):
        df3 = z[i].copy()
        df3 = df3.loc[~df3['t1'].isna()]
        df3['time'] = df3['time'] - df3['time'].min()
        df3.loc[df3['t2']>8.02, 't2'] = 8.02
        sns.kdeplot(x=df3['t2'], ax=ax, color=palette[i],
                    label=f"$j_l$: {df.iloc[i]['jl']} m/s, $j_g$: {df.iloc[i]['jg']} m/s, {df.iloc[i]['pattern_transformed']}")


    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.set_xlabel('Film thickness [mm]')
    # ax.legend(bbox_to_anchor=(1.05, 1))
    ax.legend(    loc='upper left',
    bbox_to_anchor=(1.05, 1),  # Outside the plot
    borderaxespad=0, frameon=False)

    # plt.tight_layout()
    plt.close(fig)

    return fig
