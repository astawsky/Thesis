#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations
import os


def add_lineage_regression_lines(var1, var2, fake_x):
    if ds == 'Pooled_SM':
        for dataset in pu.dataset.unique():
            ds_cond = (pu['dataset'] == dataset)
            for lin_id in pu[ds_cond].lineage_ID.unique():
                relevant = pu[ds_cond & (pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
                
                if relevant.generation.max() < 29:  # We only want long enough traces
                    continue
                
                drop_nans = relevant[[var1, var2, 'generation']].copy().dropna().sort_values('generation').reset_index(drop=True)
                
                if len(drop_nans) < 29:  # If by dropping the NaNs we do not have enough for the minimum
                    continue
                
                slope, intercept, r_value, _, _ = linregress(drop_nans[var1].values, drop_nans[var2].values)
                
                assert ~np.isnan(slope)
                
                plt.plot(fake_x, [intercept + l * slope for l in fake_x])
    else:
        for lin_id in pu.lineage_ID.unique():
            relevant = pu[(pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
            # gen_cond = (drop_nans['generation'] <= gen)
            
            if relevant.generation.max() < 29:  # We only want long enough traces
                continue
            
            drop_nans = relevant[[var1, var2, 'generation']].copy().dropna().sort_values('generation').reset_index(drop=True)
            
            if len(drop_nans) < 29:  # If by dropping the NaNs we do not have enough for the minimum
                continue
            
            slope, intercept, r_value, _, _ = linregress(drop_nans[var1].values, drop_nans[var2].values)
            
            assert ~np.isnan(slope)
            
            plt.plot(fake_x, [intercept + l * slope for l in fake_x])
            
            
def add_scatterplot_of_averages(var1, var2, pooled, df, marker='x'):
    if pooled:
        for exp, c in zip(df.experiment.unique(), cmap):
            lineages = df[(df['experiment'] == exp)].lineage_ID.unique()
            ax.scatter([df[(df['experiment'] == exp) & (df['lineage_ID'] == lin_id)][var1].mean() for lin_id in lineages[:min(len(lineages), 50)]],  # So they don't overlap too much
                       [df[(df['experiment'] == exp) & (df['lineage_ID'] == lin_id)][var2].mean() for lin_id in lineages[:min(len(lineages), 50)]],
                       marker=marker, zorder=500, label=exp.split('(')[0] if marker == 'x' else '', alpha=.5)
    else:
        ax.scatter([df[df['lineage_ID'] == lin_id][var1].mean() for lin_id in df.lineage_ID.unique()],
                   [df[df['lineage_ID'] == lin_id][var2].mean() for lin_id in df.lineage_ID.unique()], marker=marker, c='blue' if marker == 'x' else 'red', zorder=500, alpha=.5,
                   label='Trace' if marker == 'x' else 'Artificial')
            
    
def kde_scatterplot_variables(df, var1, var2, ax, line_func=None, line_label='', pooled=False, artificial=None):  # df is pu of ONE experiment
    
    sym1 = symbols['physical_units'][var1]
    sym2 = symbols['physical_units'][var2]  # r'ln$(f)$'
    
    # To speed it up we sample randomly 10,000 points
    sns.kdeplot(data=df.sample(n=10000, replace=False), x=var1, y=var2, color='gray', ax=ax)  # Do the kernel distribution approxiamtion for variables in their physical dimensions
    add_scatterplot_of_averages(var1, var2, pooled, df)  # Put the average behavior
    
    if isinstance(artificial, pd.DataFrame):
        add_scatterplot_of_averages(var1, var2, pooled, artificial, marker='^')  # How a random average behavior is supposed to act when only keeping per-cycle correlations
    
    if line_func == None:
        pass
    elif line_func == 'regression':
        x, y = df[[var1, var2]].dropna()[var1].values,df[[var1, var2]].dropna()[var2].values
        slope, intercept = linregress(x, y)[:2]
        fake_x = np.linspace(np.nanmin(x), np.nanmax(x))
        plt.plot(fake_x, intercept + slope * fake_x, color='black', ls='--', label=f'{sym2}={np.round(intercept, 2)}+{np.round(slope, 2)}*{sym1}')
    else:
        fake_x = np.linspace(np.nanmin(df[var1].values), np.nanmax(df[var1].values))
        plt.plot(fake_x, line_func(fake_x), color='black', ls='--', label=line_label)
    
    # plt.title(ds)
    ax.set_xlabel(sym1)
    ax.set_ylabel(sym2)
    
    
fig, ax = plt.subplots(tight_layout=True, figsize=[13, 6])


ds_names = cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets + sm_datasets[:-1]

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
# for ds in ds_names:
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
    
pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MC4100_37C (Tanouchi 2015)/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
art = shuffle_info(pu, True)

pu['division_ratio'] = np.log(pu['division_ratio'])
art['division_ratio'] = np.log(art['division_ratio'])

var1, var2 = 'fold_growth', 'length_birth'

num = 10
labels = np.arange(num)

data_binned, bins = pd.qcut(pu[var1].values, num, retbins=True, labels=labels)

# bin_inds = {bin: (data_binned == bin) for bin in bins}

x_binned, y_binned = [], []

for label in labels:
    indices = (data_binned == label)
    x_binned.append(pu.loc[indices][var1].mean())
    y_binned.append(pu.loc[indices][var2].mean())
    
# print(pu['div_and_fold'].mean())
# exit()

kde_scatterplot_variables(
    df=pu,
    var1=var1,
    var2=var2,
    ax=ax,
    line_func='regression',  #'regression',  #lambda x: np.log(2) / x, lambda x: -x ;;;;; None,
    pooled=False,
    artificial=art  #total_artificial
)
plt.plot(x_binned, y_binned, marker='s', label='binned', zorder=1000, alpha=.6, c='k')
plt.legend()
plt.show()
plt.close()
