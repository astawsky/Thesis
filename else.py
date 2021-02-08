#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names, shuffle_lineage_generations
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os


def plot_histograms_of_scaling_exponents(figure_folder, histogram_of_regression):
    # Histogram of H of every kind
    for kind in histogram_of_regression.kind.unique():
        print(kind)
        seaborn_preamble()
        fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
        plt.axhline(.5, ls='-', c='k', linewidth=1, zorder=1)
        
        sns.pointplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=[symbols['physical_units'][v] for v in phenotypic_variables],
                      join=False, dodge=True, capsize=.5)  # , capthick=1
        sns.boxplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=[symbols['physical_units'][v] for v in phenotypic_variables],
                    showfliers=False)
        plt.legend(title='')
        # plt.title(kind)
        plt.ylabel('')
        plt.xlabel('')
        # plt.savefig('{}/{}.png'.format(figure_folder, kind), dpi=300)
        plt.show()
        plt.close()


# def pooled_partial_corr(pu, title):
#     relevant = pu[['generationtime', 'growth_rate', 'length_birth', 'generation']].copy().dropna().sort_values('generation')
#
#     slope, intercept = linregress(relevant['generationtime'].values, relevant['growth_rate'].values)[:2]
#
#     residuals_gr = relevant['generationtime'].values - (slope * relevant['growth_rate'].values + intercept)
#
#     slope, intercept = linregress(relevant['generationtime'].values, relevant['length_birth'].values)[:2]
#
#     residuals_lb = relevant['generationtime'].values - (slope * relevant['length_birth'].values + intercept)
#
#     p = pearsonr(residuals_gr, residuals_lb)[0]
#
#     sns.regplot(residuals_gr, residuals_lb, label=f'{np.round(p, 2)}')
#     plt.title(title)
#     plt.xlabel(r'$e_{\tau, \alpha}$')
#     plt.ylabel(r'$e_{\tau, x_0}$')
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
#     plt.close()
#
#
# def lineage_specific_partial_corr(pu, title, min_length=0):
#     arr = []
#
#     for lin_id in pu['lineage_ID'].unique():
#         # print(lin_id)
#
#         relevant = pu[pu['lineage_ID'] == lin_id][['generationtime', 'growth_rate', 'length_birth', 'generation']].copy().dropna().sort_values('generation')
#
#         if len(relevant) < min_length:
#             continue
#
#         slope, intercept = linregress(relevant['generationtime'].values, relevant['growth_rate'].values)[:2]
#
#         residuals_gr = relevant['generationtime'].values - (slope * relevant['growth_rate'].values + intercept)
#
#         slope, intercept = linregress(relevant['generationtime'].values, relevant['length_birth'].values)[:2]
#
#         residuals_lb = relevant['generationtime'].values - (slope * relevant['length_birth'].values + intercept)
#
#         p = pearsonr(residuals_gr, residuals_lb)[0]
#
#         # sns.regplot(residuals_gr, residuals_lb, label=f'{np.round(p, 2)}')
#         # plt.tight_layout()
#         # plt.show()
#         # plt.close()
#
#         # print(np.round(p, 2))
#
#         arr.append(np.round(p, 2))
#
#     b = pd.DataFrame()
#     for r in arr:
#         b = b.append({r'$\rho(e_{\tau, x_0}, e_{\tau, \alpha})$': 'corr', 'partial correlation': r}, ignore_index=True)
#
#     sns.violinplot(data=b, x=r'$\rho(e_{\tau, x_0}, e_{\tau, \alpha})$', y='partial correlation')
#     plt.title(title)
#     plt.xlabel('')
#     plt.ylabel('partial_correlations')
#     plt.tight_layout()
#     plt.show()
#     plt.close()


def pooled_partial_corr(pu, title):
    
    print(pu.columns)
    relevant = pu[['generationtime', 'growth_rate', 'fold_growth', 'length_birth', 'generation']].copy().dropna().sort_values('generation')
    
    slope, intercept = linregress(relevant['generationtime'].values, relevant['fold_growth'].values)[:2]
    
    residuals_gr = relevant['generationtime'].values - (slope * relevant['fold_growth'].values + intercept)
    
    slope, intercept = linregress(relevant['fold_growth'].values, relevant['length_birth'].values)[:2]
    
    residuals_lb = relevant['fold_growth'].values - (slope * relevant['length_birth'].values + intercept)
    
    p = pearsonr(residuals_gr, residuals_lb)[0]
    
    sns.regplot(residuals_gr, residuals_lb, label=f'{np.round(p, 2)}')
    plt.title(title)
    plt.xlabel(r'$e_{\alpha, \tau}$')
    plt.ylabel(r'$e_{\alpha, x_0}$')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()


def lineage_specific_partial_corr(pu, title, min_length=0):
    arr = []
    
    for lin_id in pu['lineage_ID'].unique():
        # print(lin_id)
        
        relevant = pu[pu['lineage_ID'] == lin_id][['generationtime', 'growth_rate', 'fold_growth', 'length_birth', 'generation']].copy().dropna().sort_values('generation')
        
        if len(relevant) < min_length:
            continue
        
        slope, intercept = linregress(relevant['generationtime'].values, relevant['fold_growth'].values)[:2]
        
        residuals_gr = relevant['generationtime'].values - (slope * relevant['fold_growth'].values + intercept)
        
        slope, intercept = linregress(relevant['fold_growth'].values, relevant['length_birth'].values)[:2]
        
        residuals_lb = relevant['fold_growth'].values - (slope * relevant['length_birth'].values + intercept)
        
        p = pearsonr(residuals_gr, residuals_lb)[0]
        
        # sns.regplot(residuals_gr, residuals_lb, label=f'{np.round(p, 2)}')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # print(np.round(p, 2))
        
        arr.append(np.round(p, 2))
    
    b = pd.DataFrame()
    for r in arr:
        b = b.append(pd.DataFrame({r'$\rho(e_{\phi, x_0}, e_{\phi, \tau})$': 'corr', 'partial correlation': r}, index=[len(b)]), ignore_index=True)
        
    # print(b)
    # exit()
    
    sns.violinplot(x=b[list(b.columns)[0]], y=b['partial correlation'])
    plt.title(title)
    plt.xlabel('')
    plt.ylabel('partial_correlations')
    plt.tight_layout()
    plt.show()
    plt.close()


def check_simulation(samples=1000, min_length=0):
    variables = np.random.multivariate_normal([0, 0, 0],
                                              [[1, 0, -.5],
                                               [0, 1, -.3],
                                               [-.5, -.3, 1]], samples)
    
    x = variables[:, 0]
    a = variables[:, 1]
    t = variables[:, 2]
    
    # This is to show that via a multivariate Gaussian there is no way to get as high a correlation for fold_growth
    # and size than generationtime and size if generationtime and growth rate are negatively correlated
    for var1, var2, label in [(x, a, 'x, a'), (x, t, 'x, t'), (a, t, 'a, t'), (x, a * t, 'x, a*t')]:
        sns.regplot(var1, var2, label=f'{label}: {np.round(pearsonr(var1, var2)[0], 2)}')
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()
    
    pu = pd.DataFrame.from_dict({num: {
        'generationtime': t[num],
        'growth_rate': a[num],
        'fold_growth': a[num]*t[num],
        'length_birth': x[num],
        'generation': num
    } for num in np.arange(samples)}, 'index')

    pooled_partial_corr(pu, 'Simulation')
    # lineage_specific_partial_corr(pu, 'Simulation', min_length)


# for data_origin in ['Pooled_SM']:
#     print(data_origin)
#
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#
#     # create_folder(current_dir + '/Dataframes')
#     # create_folder(current_dir + '/Dataframes/' + data_origin)
#
#     # create_folder(current_dir + '/Scaling Exponents')
#     # create_folder(current_dir + '/LogLog Recreation')
#     # create_folder(current_dir + '/Scaling Exponents/' + data_origin)
#     # create_folder(current_dir + '/LogLog Recreation/' + data_origin)
#
#     processed_data = os.path.dirname(current_dir) + '/Datasets/' + data_origin + '/ProcessedData/'
#
#     """
#     data_origin ==> Name of the dataset we are analysing
#     raw_data ==> Where the folder containing the raw data for this dataset is
#     processed_data ==> The folder we will put the processed data in
#     """
#
#     args = {
#         'data_origin': data_origin,
#         # Data singularities, long traces with significant filamentation, sudden drop-offs
#         'Scaling_Exponents': current_dir + '/Scaling Exponents/' + data_origin,
#         'LogLog_Recreation': current_dir + '/LogLog Recreation/' + data_origin,
#         # 'dataframes': current_dir + '/Dataframes/' + data_origin,
#         'pu': processed_data + '/physical_units.csv'
#     }
#
#     # main(args)
#
#     df = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/DFA/Scaling Exponents/Pooled_SM/scaling_exponents.csv')  # '{}/scaling_exponents.csv'.format(args['Scaling_Exponents'])
#
#     plot_histograms_of_scaling_exponents(args['Scaling_Exponents'], df)

for ds in ['Pooled_SM']:
    print(ds)
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/physical_units.csv')
    # ta = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()

    table = pd.DataFrame()
    meaned = pd.DataFrame()
    
    for lin_id in pu.lineage_ID.unique():
        relevant = pu[pu['lineage_ID'] == lin_id].copy()
        
        if len(relevant) < 30:  # Need ones that are long enough
            continue
            
        for vari in phenotypic_variables:
            table = table.append({
                'variable': vari,
                'std': relevant[vari].std(),
                'mean': relevant[vari].mean(),
                'cv': relevant[vari].std() / relevant[vari].mean(),
                'lineage_ID': lin_id
            }, ignore_index=True)

    for vari in phenotypic_variables:
        meaned = meaned.append({
            'variable': vari,
            'cv': table[table['variable'] == vari]['cv'].mean()
        }, ignore_index=True)
        
    print(meaned)
    
    # sns.violinplot(data=table, x='variable', y='cv')
    # plt.show()
    # plt.close()

exit()

# check_simulation(samples=1000, min_length=0)

for ds in dataset_names:
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/physical_units.csv')
    
    pooled_partial_corr(pu, ds)
    lineage_specific_partial_corr(pu, ds, min_length=20)
    pooled_partial_corr(pu, ds)

exit()

# pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MC4100_25C (Tanouchi 2015)/ProcessedData/physical_units.csv')

relevant = pu[['generationtime', 'growth_rate', 'length_birth', 'generation']].copy().dropna().sort_values('generation')

slope, intercept = linregress(relevant['generationtime'].values, relevant['growth_rate'].values)[:2]

residuals_gr = relevant['generationtime'].values - (slope * relevant['growth_rate'].values + intercept)

slope, intercept = linregress(relevant['generationtime'].values, relevant['length_birth'].values)[:2]

residuals_lb = relevant['generationtime'].values - (slope * relevant['length_birth'].values + intercept)

p = pearsonr(residuals_gr, residuals_lb)[0]

sns.regplot(residuals_gr, residuals_lb, label=f'{np.round(p, 2)}')
plt.legend()
plt.tight_layout()
plt.show()
plt.close()
#
# print(np.round(p, 2))
#
# exit()
#
#
variables = np.random.multivariate_normal([0, 0, 0],
                                          [[1, 0, -.5],
                                           [0, 1, -.3],
                                           [-.5, -.3, 1]], 1000)

x = variables[:, 0]
a = variables[:, 1]
t = variables[:, 2]

# This is to show that via a multivariate Gaussian there is no way to get as high a correlation for fold_growth
# and size than generationtime and size if generationtime and growth rate are negatively correlated
for var1, var2, label in [(x, a, 'x, a'), (x, t, 'x, t'), (a, t, 'a, t'), (x, a * t, 'x, a*t')]:
    sns.regplot(var1, var2, label=f'{label}: {np.round(pearsonr(var1, var2)[0], 2)}')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()

exit()
