#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names
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
        first = [df[df['lineage_ID'] == lin_id][var1].mean() for lin_id in df.lineage_ID.unique()]
        second = [df[df['lineage_ID'] == lin_id][var2].mean() for lin_id in df.lineage_ID.unique()]
        ax.scatter(first, second, marker=marker, c='blue' if marker == 'x' else 'red', zorder=500, alpha=.5 if marker == 'x' else .2,
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
            
    
def kde_scatterplot_variables(df, var1, var2, num, ax, line_func=None, line_label='', pooled=False, artificial=None):  # df is pu of ONE experiment
    
    sym1 = symbols['physical_units'][var1]
    sym2 = symbols['physical_units'][var2]  # r'ln$(f)$'
    
    # To speed it up we sample randomly 10,000 points
    sns.kdeplot(data=df.sample(n=10000, replace=False), x=var1, y=var2, color='gray', ax=ax)  # Do the kernel distribution approxiamtion for variables in their physical dimensions
    add_scatterplot_of_averages(var1, var2, pooled, df, ax, marker='x')  # Put the average behavior
    
    if isinstance(artificial, pd.DataFrame):
        add_scatterplot_of_averages(var1, var2, pooled, artificial, ax, marker='^')  # How a random average behavior is supposed to act when only keeping per-cycle correlations
    
    if line_func == None:
        pass
    elif line_func == 'regression':
        x, y = df[[var1, var2]].dropna()[var1].values,df[[var1, var2]].dropna()[var2].values
        slope, intercept = linregress(x, y)[:2]
        fake_x = np.linspace(np.nanmin(x), np.nanmax(x))
        ax.plot(fake_x, intercept + slope * fake_x, color='black', ls='--', label=f'{sym2}={np.round(intercept, 2)}+{np.round(slope, 2)}*{sym1}')
    else:
        fake_x = np.linspace(np.nanmin(df[var1].values), np.nanmax(df[var1].values))
        ax.plot(fake_x, line_func(fake_x), color='black', ls='--', label=line_label)
    
    # plt.title(ds)
    plot_binned_data(df, var1, var2, num, ax)
    ax.set_xlabel(sym1)
    ax.set_ylabel(sym2)
    # ax.legend(title='')

    no_nans = df[[var1, var2]].copy().dropna()

    no_nans_art = artificial[[var1, var2]].copy().dropna()

    print('kde scatter plot')
    print(f'pooled correlation: {pearsonr(no_nans[var1].values, no_nans[var2].values)[0]} = {pearsonr(no_nans_art[var1].values, no_nans_art[var2].values)[0]}')
    print('-'*200)
    
    
def plot_pair_scatterplots(df, var1, var2, ax):
    x_a = [df[(df['trap_ID'] == trap_id) & (df['trace'] == 'A')][var1].mean() for trap_id in df.trap_ID.unique()]
    y_b = [df[(df['trap_ID'] == trap_id) & (df['trace'] == 'B')][var2].mean() for trap_id in df.trap_ID.unique()]
    
    ax.grid(True)
    ax.scatter(x_a, y_b)
    ax.set_xlabel(symbols['time_averages'][var1] if var1 != 'division_ratio' else r'$\ln(\overline{f}^A)$')
    ax.set_ylabel(symbols['time_averages'][var2] if var2 != 'division_ratio' else r'$\ln(\overline{f}^B)$')
    
    print(f'pair scatterplot {var1} {var2}: {pearsonr(x_a, y_b)[0]}')


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

# total_ta = pd.DataFrame()
# total = pd.DataFrame()
# total_artificial = pd.DataFrame()
# for ds in dataset_names[:-1]:
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
    
pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
pu['division_ratio'] = np.log(pu['division_ratio'])
art = shuffle_info(pu, False)

num = 10
    
# print(pu['div_and_fold'].mean())
# exit()

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

axes[0, 0].set_xlim([.76, 2])
axes[0, 0].set_ylim([.19, .945])

axes[0, 1].set_xlim([.364, .725])
axes[0, 1].set_ylim([.398, .692])

axes[0, 2].set_xlim([.93, 1.646])
axes[0, 2].set_ylim([.922, 1.597])

###################################

axes[1, 0].set_xlim([-.798, -.39])
axes[1, 0].set_ylim([.21, 1.078])

axes[1, 1].set_xlim([.5453, .742])
axes[1, 1].set_ylim([.539, .742])

axes[1, 2].set_xlim([-.6802, -.4895])
axes[1, 2].set_ylim([-.6614, -.5083])

kde_scatterplot_variables(
    df=pu,
    var1='growth_rate',
    var2='generationtime',
    num=num,
    ax=axes[0, 0],
    line_func=None,  #'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
    pooled=False,
    artificial=art
)

kde_scatterplot_variables(
    df=pu,
    var1='division_ratio',
    var2='fold_growth',
    num=num,
    ax=axes[1, 0],
    line_func=None,  #'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
    pooled=False,
    artificial=art
)

plot_pair_scatterplots(pu, 'generationtime', 'generationtime', axes[0, 1])
plot_pair_scatterplots(pu, 'growth_rate', 'growth_rate', axes[0, 2])
plot_pair_scatterplots(pu, 'fold_growth', 'fold_growth', axes[1, 1])
plot_pair_scatterplots(pu, 'division_ratio', 'division_ratio', axes[1, 2])
# plt.legend()
# plt.savefig('Figures/lineage information and micro environment with titles', dpi=300)
plt.show()
plt.close()
