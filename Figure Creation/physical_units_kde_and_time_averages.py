#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names, pool_experiments
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os
            
            
def add_scatterplot_of_averages(var1, var2, pooled, df, ax, marker='x'):
    if pooled:
        for exp, c in zip(df.experiment.unique(), cmap):
            lineages = df[(df['experiment'] == exp)].lineage_ID.unique()
            ax.scatter([df[(df['experiment'] == exp) & (df['lineage_ID'] == lin_id)][var1].mean() for lin_id in lineages[:min(len(lineages), 50)]],  # So they don't overlap too much
                       [df[(df['experiment'] == exp) & (df['lineage_ID'] == lin_id)][var2].mean() for lin_id in lineages[:min(len(lineages), 50)]],
                       marker=marker, zorder=500, label=exp.split('(')[0] if marker == 'x' else '', alpha=.5)
    else:
        first = [df[df['lineage_ID'] == lin_id][var1].mean() for lin_id in df.lineage_ID.unique() if len(df[df['lineage_ID'] == lin_id]) > 6]
        second = [df[df['lineage_ID'] == lin_id][var2].mean() for lin_id in df.lineage_ID.unique() if len(df[df['lineage_ID'] == lin_id]) > 6]
        ax.scatter(first, second, marker=marker, c=cmap[0] if marker == 'x' else cmap[1], zorder=500, alpha=.5 if marker == 'x' else .2,
                   label='Trace' if marker == 'x' else 'Artificial')
        print(f'Trace averages correlation: {pearsonr(first, second)[0]}' if marker == 'x' else f'Artificial averages correlation: {pearsonr(first, second)[0]}')
        
    
def plot_binned_data(df, var1, var2, num, ax):
    labels = np.arange(num)
    data_binned, bins = pd.qcut(df[var1].values, num, retbins=True, labels=labels)
    
    # bin_inds = {bin: (data_binned == bin) for bin in bins}
    
    x_binned, y_binned = [], []
    
    for label in labels:
        indices = (data_binned == label)
        x_binned.append(df.loc[indices][var1].mean())
        y_binned.append(df.loc[indices][var2].mean())
        
    ax.plot(x_binned, y_binned, marker='s', label='binned', zorder=1000, alpha=.6, c='k')
            
    
def kde_scatterplot_variables(df, var1, var2, num, ax, line_func=[], line_label='', pooled=False, artificial=None, sym1=None, sym2=None):  # df is pu of ONE experiment
    
    if sym1 == None:
        sym1 = symbols['physical_units'][var1]# if var1 != 'division_ratio' else r'$\ln(f)$'
    if sym2 == None:
        sym2 = symbols['physical_units'][var2]
    
    # To speed it up we sample randomly 1,000 points
    sns.kdeplot(data=df.sample(n=1000, replace=False), x=var1, y=var2, color='gray', ax=ax)  # Do the kernel distribution approxiamtion for variables in their physical dimensions
    add_scatterplot_of_averages(var1, var2, pooled, df, ax, marker='x')  # Put the average behavior
    
    if isinstance(artificial, pd.DataFrame):
        add_scatterplot_of_averages(var1, var2, pooled, artificial, ax, marker='^')  # How a random average behavior is supposed to act when only keeping per-cycle correlations

        no_nans_art = artificial[[var1, var2]].copy().dropna()
        print(f'pooled correlation (Artificial): {pearsonr(no_nans_art[var1].values, no_nans_art[var2].values)[0]}')
    
    if len(line_func) == 0:
        pass
    elif line_func == 'regression':
        x, y = df[[var1, var2]].dropna()[var1].values,df[[var1, var2]].dropna()[var2].values
        slope, intercept = linregress(x, y)[:2]
        fake_x = np.linspace(np.nanmin(x), np.nanmax(x))
        ax.plot(fake_x, intercept + slope * fake_x, color='black', ls='--', label=f'{sym2}={np.round(intercept, 2)}+{np.round(slope, 2)}*{sym1}')
    else:
        for count, func in enumerate(line_func):
            print(count)
            print(func)
            fake_x = np.linspace(np.nanmin(df[var1].values), np.nanmax(df[var1].values))
            ax.plot(fake_x, func(fake_x), color='black' if count == 0 else 'gray', ls='--', label=line_label)
    
    # plt.title(ds)
    plot_binned_data(df, var1, var2, num, ax)
    ax.set_xlabel(sym1)
    ax.set_ylabel(sym2)
    # ax.legend(title='')

    no_nans = df[[var1, var2]].copy().dropna()

    print('kde scatter plot')
    print(f'pooled correlation (Trace): {pearsonr(no_nans[var1].values, no_nans[var2].values)[0]}')
    print('-'*200)
    
    
def plot_pair_scatterplots(df, var1, var2, ax, sym1=None, sym2=None):
    # x_a = [df[(df['trap_ID'] == trap_id) & (df['trace'] == 'A') & (df['dataset'] == 'NL')][var1].mean() for trap_id in df[(df['dataset'] == 'NL')].trap_ID.unique()]
    # y_b = [df[(df['trap_ID'] == trap_id) & (df['trace'] == 'B') & (df['dataset'] == 'NL')][var2].mean() for trap_id in df[(df['dataset'] == 'NL')].trap_ID.unique()]
    x_a = []
    y_b = []
    
    # diagonal_length = []
    # offdiagonal_length = []

    for trap_id in df[(df['dataset'] == 'NL')].trap_ID.unique():
        lin_a = df[(df['trap_ID'] == trap_id) & (df['trace'] == 'A') & (df['dataset'] == 'NL')].copy()
        lin_b = df[(df['trap_ID'] == trap_id) & (df['trace'] == 'B') & (df['dataset'] == 'NL')].copy()
        
        if (len(lin_a) < 7) or (len(lin_b) < 7):
            continue
        x_a.append(lin_a[var1].mean())
        y_b.append(lin_b[var2].mean())
    
    x_a = np.array(x_a)
    y_b = np.array(y_b)
    
    ax.grid(True)
    ax.scatter(x_a, y_b)
    
    diag_spead = np.std((x_a + y_b) / 2)
    off_spread = np.std((np.append((x_a - np.sqrt(x_a * y_b)), (y_b - np.sqrt(x_a * y_b)))))
    center = np.mean(np.append(x_a, y_b))
    
    ax.plot(np.linspace(center - diag_spead, center + diag_spead), np.linspace(center - diag_spead, center + diag_spead), ls='-', c='k', linewidth=3)
    ax.plot(np.linspace(center - off_spread, center + off_spread), - np.linspace(center - off_spread, center + off_spread) + 2 * center, ls='-', c='k', linewidth=3)
    
    # slope, intercept = linregress(x_a, y_b)[:2]
    # low_x, low_y = np.mean(x_a) - np.std((np.array(y_b)+np.array(x_a))/2), np.mean(y_b) - np.std(np.array(y_b)-(intercept+slope*np.array(x_a)))
    # high_x, high_y = np.mean(x_a) + np.std((np.array(y_b)+np.array(x_a))/2), np.mean(y_b) + np.std(np.array(y_b)-(intercept+slope*np.array(x_a)))
    #
    # diag_low, diag_high = np.mean(x_a) - np.std((np.array(y_b)+np.array(x_a))/2), np.mean(x_a) + np.std((np.array(y_b)+np.array(x_a))/2)
    #
    # new_int = (np.mean(y_b) + np.mean(x_a) / slope)
    #
    # ax.plot(np.linspace(low_x, high_x), intercept + slope * np.linspace(low_x, high_x), ls='-', c='k', linewidth=3)
    #
    # ax.plot(np.linspace((-low_y + new_int)*slope, (-high_y + new_int)*slope), (np.mean(y_b) + np.mean(x_a) / slope) - np.linspace((-low_y + new_int)*slope, (-high_y + new_int)*slope) / slope, ls='-', c='k', linewidth=3)
    
    x_a, y_b = list(x_a), list(y_b)
    
    minimum, maximum = np.min(x_a + y_b) - (.2*np.std(x_a + y_b)), np.max(x_a + y_b) + (.2*np.std(x_a + y_b))
    ax.set_xlim([minimum, maximum])
    ax.set_ylim([minimum, maximum])
    ax.set_xticks(np.round(np.linspace(minimum, maximum, 4), 2))
    ax.set_yticks(np.round(np.linspace(minimum, maximum, 4), 2))
    if sym1 == None:
        sym1 = symbols['time_averages'][var1] + r'$^{\, A}$'
    if sym2 == None:
        sym2 = symbols['time_averages'][var2] + r'$^{\, B}$'
    ax.set_xlabel(sym1)
    ax.set_ylabel(sym2)
    
    print(f'pair scatterplot {var1} {var2}: {pearsonr(x_a, y_b)[0]}')
    
    
#########################################
# What appears in the main text #
#########################################
    
    
def alpha_tau_script():
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['fold_growth'] = np.log(pu['fold_growth'])
    art = shuffle_info(pu, False)
    
    num = 10

    scale = 1.5

    sns.set_context('paper', font_scale=1 * scale)
    sns.set_style("ticks", {'axes.grid': False})

    fig, axes = plt.subplots(2, 3, tight_layout=True, figsize=[6.5 * scale, 4.2 * scale])

    axes[0, 0].set_title('A', x=-.2, fontsize='xx-large')
    axes[0, 1].set_title('B', x=-.2, fontsize='xx-large')
    axes[0, 2].set_title('C', x=-.2, fontsize='xx-large')
    axes[1, 0].set_title('D', x=-.2, fontsize='xx-large')
    axes[1, 1].set_title('E', x=-.2, fontsize='xx-large')
    axes[1, 2].set_title('F', x=-.2, fontsize='xx-large')

    axes[0, 2].set_xlim([.76, 2])
    axes[0, 2].set_ylim([.19, 1])
    axes[1, 2].set_xlim([.4, .63])
    axes[1, 2].set_ylim([np.exp(.3), np.exp(1.15)])

    kde_scatterplot_variables(
        df=pu,
        var1='growth_rate',
        var2='generationtime',
        num=num,
        ax=axes[0, 2],
        line_func=[lambda x: np.log(2) / x, lambda x: pu['fold_growth'].mean() / x ],  # 'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
        pooled=False,
        artificial=art
    )

    changed = pu.copy()
    changed['fold_growth'] = np.exp(changed['fold_growth'])

    kde_scatterplot_variables(
        df=changed,
        var1='division_ratio',
        var2='fold_growth',
        num=num,
        ax=axes[1, 2],
        line_func=[lambda x: 1 / x, lambda x: np.mean(np.exp(pu['fold_growth']) * np.mean(pu['division_ratio'])) / x],  # 'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
        pooled=False,
        artificial=shuffle_info(changed, False),
        sym2=r'$e^{\phi}$'
    )

    plot_pair_scatterplots(pu, 'generationtime', 'generationtime', axes[0, 1])
    plot_pair_scatterplots(pu, 'growth_rate', 'growth_rate', axes[0, 0])
    plot_pair_scatterplots(pu, 'fold_growth', 'fold_growth', axes[1, 0])
    plot_pair_scatterplots(pu, 'division_ratio', 'division_ratio', axes[1, 1])  # , sym1=r'$\overline{\ln(f)}^{\, A}$', sym2=r'$\overline{\ln(f)}^{\, B}$')
    # plt.legend()
    # plt.savefig('Figures/lineage information and micro environment with titles', dpi=300)
    plt.show()
    plt.close()
    
    
def size_script():
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    pu['fold_growth'] = np.exp(pu['fold_growth'])
    art = shuffle_info(pu, False)
    
    num = 10
    
    scale = 1.5
    
    sns.set_context('paper', font_scale=1 * scale)
    sns.set_style("ticks", {'axes.grid': False})
    
    fig, axes = plt.subplots(1, 3, tight_layout=True, figsize=[6.5 * scale, 2.1 * scale])
    
    axes[0].set_title('A', x=-.2, fontsize='xx-large')
    axes[1].set_title('B', x=-.2, fontsize='xx-large')
    axes[2].set_title('C', x=-.2, fontsize='xx-large')
    
    axes[0].set_xlim([1.427, 4.278])
    axes[0].set_ylim([2.8, 7.9])

    axes[1].set_xlim([1.404, 4.221])
    axes[1].set_ylim([.78, 4.64])

    axes[2].set_xlim([1.537, 4.328])
    axes[2].set_ylim([1.34, 3.088])
    
    kde_scatterplot_variables(
        df=pu,
        var1='length_birth',
        var2='length_final',
        num=num,
        ax=axes[0],
        line_func=[lambda x: 2 * x],  # 'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
        pooled=False,
        artificial=art
    )
    
    print(pu['length_birth'].values - (np.mean(pu['length_birth'].values) * (-1 + np.exp(pu['fold_growth'].mean())))/2)
    
    kde_scatterplot_variables(
        df=pu,
        var1='length_birth',
        var2='added_length',
        num=num,
        ax=axes[1],
        line_func=[lambda x: x,
                   lambda x: x - (np.nanmean(x) * (-1 + np.exp(pu['fold_growth'].mean())))/2],  # 'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
        pooled=False,
        artificial=art
    )
    
    kde_scatterplot_variables(
        df=pu,
        var1='length_birth',
        var2='fold_growth',
        num=num,
        ax=axes[2],
        line_func=[lambda x: 2 * np.array([1 for _ in np.arange(len(x))]),
                   lambda x: np.nanmean(pu['fold_growth'].values) * np.array([1 for _ in np.arange(len(x))])],  # 'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
        pooled=False,
        artificial=art,
        sym2=r'$e^{\phi}$'
    )
    
    # plt.legend()
    plt.savefig('Figures/size lineage information and micro environment with titles', dpi=300)
    # plt.show()
    plt.close()
    
    
#########################################
# What appears in the appendix #
#########################################
    
    
def size_and_rest_script():
    
    """
    Here we show that average growth rate or inter-div. time and size are uncorrelated between lineage averages.
    
    :return: plots
    """

    # pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['fold_growth'] = np.exp(pu['fold_growth'])
    # art = shuffle_info(pu, False)

    total = pd.DataFrame()
    total_artificial = pd.DataFrame()

    for name, group in zip(['Wang 2010, CGSC 6300', 'Wang 2010, lexA3', 'Susman and Kohram', 'Vashistha 2021', 'Tanouchi 25C', 'Tanouchi 27C', 'Tanouchi 37C'],
                           [cgsc_6300_wang_exps, lexA3_wang_exps, mm_datasets, ['Pooled_SM'], ['MC4100_25C (Tanouchi 2015)'], ['MC4100_27C (Tanouchi 2015)'], ['MC4100_37C (Tanouchi 2015)']]):
        print(name, group)
    
        # df = pool_experiments(group, name, dimensions='pu')
    
        df = pool_experiments(group, name, dimensions='ta')[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
        df = df[df['max_gen'] > 7]
    
        art = shuffle_info(df, False if name == 'Pooled_SM' else True)
    
        df['experiment'] = name
        art['experiment'] = name
        df['division_ratio'] = df['division_ratio']
        art['division_ratio'] = art['division_ratio']
    
        total = total.append(df, ignore_index=True)
        total_artificial = total_artificial.append(art, ignore_index=True)

    num = 10

    scale = 1.5

    sns.set_context('paper', font_scale=1 * scale)
    sns.set_style("ticks", {'axes.grid': False})

    # fig, ax = plt.subplots(tight_layout=True, figsize=[6.5, 6.5])
    fig, axes = plt.subplots(2, 2, tight_layout=True, figsize=[6.5, 7.5])

    axes[0, 0].set_title('A', x=-.2, fontsize='xx-large')
    axes[0, 1].set_title('B', x=-.2, fontsize='xx-large')
    axes[1, 0].set_title('C', x=-.2, fontsize='xx-large')
    axes[1, 1].set_title('D', x=-.2, fontsize='xx-large')

    kde_scatterplot_variables(
        df=total,
        var1='growth_rate',
        var2='length_birth',
        num=num,
        ax=axes[0, 0],
        line_func=[],
        pooled=True
    )

    kde_scatterplot_variables(
        df=total,
        var1='generationtime',
        var2='length_birth',
        num=num,
        ax=axes[0, 1],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='division_ratio',
        var2='length_birth',
        num=num,
        ax=axes[1, 0],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='div_and_fold',
        var2='length_birth',
        num=num,
        ax=axes[1, 1],
        line_func=[],
        pooled=True
    )

    axes[0, 1].legend(title='', fontsize='xx-small', markerscale=1, loc='upper right') # , bbox_to_anchor=(4, 2)
    # plt.savefig('Figures/rest of size correlations with experiments', dpi=300)
    plt.show()
    plt.close()
    
    
def manifold_across_experiments():
    """
        Here we show that average growth rate or inter-div. time and size are uncorrelated between lineage averages.

        :return: plots
        """
    
    # pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['fold_growth'] = np.exp(pu['fold_growth'])
    # art = shuffle_info(pu, False)
    
    total = pd.DataFrame()
    total_artificial = pd.DataFrame()
    
    for name, group in zip(['Wang 2010, CGSC 6300', 'Wang 2010, lexA3', 'Susman and Kohram', 'Vashistha 2021', 'Tanouchi 25C', 'Tanouchi 27C', 'Tanouchi 37C'],
                           [cgsc_6300_wang_exps, lexA3_wang_exps, mm_datasets, ['Pooled_SM'], ['MC4100_25C (Tanouchi 2015)'], ['MC4100_27C (Tanouchi 2015)'], ['MC4100_37C (Tanouchi 2015)']]):
        print(name, group)
        
        # df = pool_experiments(group, name, dimensions='pu')
        
        df = pool_experiments(group, name, dimensions='ta')[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
        df = df[df['max_gen'] > 7]
        
        art = shuffle_info(df, False if name == 'Pooled_SM' else True)
        
        df['experiment'] = name
        art['experiment'] = name
        df['division_ratio'] = df['division_ratio']
        art['division_ratio'] = art['division_ratio']
        
        total = total.append(df, ignore_index=True)
        total_artificial = total_artificial.append(art, ignore_index=True)
    
    num = 10
    
    scale = 1.5
    
    sns.set_context('paper', font_scale=1 * scale)
    sns.set_style("ticks", {'axes.grid': False})
    
    # fig, ax = plt.subplots(tight_layout=True, figsize=[6.5, 6.5])
    fig, axes = plt.subplots(2, 2, tight_layout=True, figsize=[6.5, 7.5])
    
    axes[0, 0].set_title('A', x=-.2, fontsize='xx-large')
    axes[0, 1].set_title('B', x=-.2, fontsize='xx-large')
    axes[1, 0].set_title('C', x=-.2, fontsize='xx-large')
    axes[1, 1].set_title('D', x=-.2, fontsize='xx-large')
    
    kde_scatterplot_variables(
        df=total,
        var1='growth_rate',
        var2='div_and_fold',
        num=num,
        ax=axes[0, 0],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='generationtime',
        var2='div_and_fold',
        num=num,
        ax=axes[0, 1],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='division_ratio',
        var2='div_and_fold',
        num=num,
        ax=axes[1, 0],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='fold_growth',
        var2='div_and_fold',
        num=num,
        ax=axes[1, 1],
        line_func=[],
        pooled=True
    )
    
    axes[0, 1].legend(title='', fontsize='xx-small', markerscale=1, loc='upper right')  # , bbox_to_anchor=(4, 2)
    # plt.savefig('Figures/rest of size correlations with experiments', dpi=300)
    plt.show()
    plt.close()
    
    
def between_experiments_correlation():
    """
        Here we show that average growth rate or inter-div. time and size are uncorrelated between lineage averages.

        :return: plots
        """
    
    # pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['fold_growth'] = np.exp(pu['fold_growth'])
    # art = shuffle_info(pu, False)
    
    total = pd.DataFrame()
    total_artificial = pd.DataFrame()
    
    for name, group in zip(['Wang 2010, CGSC 6300', 'Wang 2010, lexA3', 'Susman and Kohram', 'Vashistha 2021', 'Tanouchi 25C', 'Tanouchi 27C', 'Tanouchi 37C'],
                           [cgsc_6300_wang_exps, lexA3_wang_exps, mm_datasets, ['Pooled_SM'], ['MC4100_25C (Tanouchi 2015)'], ['MC4100_27C (Tanouchi 2015)'], ['MC4100_37C (Tanouchi 2015)']]):
        print(name, group)
        
        # df = pool_experiments(group, name, dimensions='pu')
        
        df = pool_experiments(group, name, dimensions='ta')[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
        df = df[df['max_gen'] > 7]
        
        art = shuffle_info(df, False if name == 'Pooled_SM' else True)
        
        df['experiment'] = name
        art['experiment'] = name
        
        total = total.append(df, ignore_index=True)
        total_artificial = total_artificial.append(art, ignore_index=True)
    
    num = 10
    
    scale = 1.5
    
    sns.set_context('paper', font_scale=1 * scale)
    sns.set_style("ticks", {'axes.grid': False})
    
    # fig, ax = plt.subplots(tight_layout=True, figsize=[6.5, 6.5])
    fig, axes = plt.subplots(2, 2, tight_layout=True, figsize=[6.5, 7.5])
    
    axes[0, 0].set_title('A', x=-.2, fontsize='xx-large')
    axes[0, 1].set_title('B', x=-.2, fontsize='xx-large')
    axes[1, 0].set_title('C', x=-.2, fontsize='xx-large')
    axes[1, 1].set_title('D', x=-.2, fontsize='xx-large')
    
    kde_scatterplot_variables(
        df=total,
        var1='growth_rate',
        var2='division_ratio',
        num=num,
        ax=axes[0, 0],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='generationtime',
        var2='division_ratio',
        num=num,
        ax=axes[0, 1],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='length_birth',
        var2='added_length',
        num=num,
        ax=axes[1, 0],
        line_func=[],
        pooled=True
    )
    
    kde_scatterplot_variables(
        df=total,
        var1='division_ratio',
        var2='added_length',
        num=num,
        ax=axes[1, 1],
        line_func=[],
        pooled=True
    )
    
    axes[0, 1].legend(title='', fontsize='xx-small', markerscale=1, loc='upper right')  # , bbox_to_anchor=(4, 2)
    # plt.savefig('Figures/rest of size correlations with experiments', dpi=300)
    plt.show()
    plt.close()


between_experiments_correlation()
# manifold_across_experiments()
# size_and_rest_script()
# alpha_tau_script()
# size_script()
exit()


# ds_names = cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets + sm_datasets[:-1]
#
# exit()
#
#
# for ds in ds_names:
#     print(ds)
#     kde_scatterplot_variables(
#         pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'),
#         var1='generationtime',
#         var2='growth_rate',
#         ax=ax,
#         line_func=None,
#         pooled=False
#     )
#
#     plt.show()
#     plt.close()
#
# exit()

# for ds in dataset_names:
#     if ds == 'Pooled_SM':
#         continue
#
#     print(ds)
#
#     pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
#     art = shuffle_info(pu, False if ds in sm_datasets else True)
#
#     # if ds not in ['MG1655_inLB_LongTraces'] + sm_datasets[:-1]:
#     #     ta = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
#     # else:
#     #     print('TA not available before')
#     #     ta = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
#     #     ta.to_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
#
#     # ta = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
#     # ta.to_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
#     # ta['experiment'] = ds
#     pu['experiment'] = ds
#     art['experiment'] = ds
#     pu['division_ratio'] = np.log(pu['division_ratio'])
#     art['division_ratio'] = np.log(art['division_ratio'])
#
#     total = total.append(pu, ignore_index=True)
#     total_artificial = total_artificial.append(art, ignore_index=True)
#
#     # total_ta = total_ta.append(ta, ignore_index=True)
    

