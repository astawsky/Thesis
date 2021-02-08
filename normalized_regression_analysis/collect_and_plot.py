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

ds_names = [sm_datasets[-1]] + cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets
ds_names.remove('Maryam_LongTraces')

same_cell_variable_combos = list(combinations(['length_birth', 'generationtime', 'growth_rate', 'division_ratio'], 2))
same_cell_variable_combos = same_cell_variable_combos + [('fold_growth', 'length_birth'), ('length_final', 'length_birth'), ('fold_growth', 'division_ratio'), ('added_length', 'division_ratio')]

# pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MC4100_25C (Tanouchi 2015)/ProcessedData/physical_units.csv')
# print(pu['length_birth'])
# pu['length_birth'] = np.log(pu['length_birth'])
# print(pu['length_birth'])

# print(linregress(pu.fold_growth.values, pu.length_birth.values))
#
# sns.regplot(data=pu, x='fold_growth', y='length_birth')
# plt.show()
# plt.close()
# exit()
#
# regression_dataframe = pd.DataFrame()
#
#
# def get_regression_data():
#     for n in np.arange(30, len(relevant) + 2, 5):
#         for var1, var2 in same_cell_variable_combos:
#
#             drop_nans = relevant[[var1, var2, 'generation']].copy().dropna().sort_values('generation').reset_index(drop=True)
#             # gen_cond = (drop_nans['generation'] <= gen)
#
#             if len(drop_nans) < n:  # If by dropping the NaNs we do not have enough for the minimum
#                 continue
#
#             x, y = drop_nans[var1].iloc[:n], drop_nans[var2].iloc[:n]
#             x, y = ((x - x.mean()) / x.std()).values, ((y - y.mean()) / y.std()).values
#
#             slope, intercept, r_value, _, _ = linregress(x, y)
#
#             assert ~np.isnan(slope)
#             assert ~np.isnan(intercept)
#             assert ~np.isnan(r_value)
#             assert (len(x) % 5 == 0)
#             assert (n % 5 == 0)
#
#             temp_dict.update({
#                 len(temp_dict): {
#                     'slope': slope,
#                     'intercept': intercept,
#                     'r_value': r_value,
#                     'lin_id': lin_id,
#                     'ds': ds,
#                     'n': n,
#                     'var1': var1,
#                     'var2': var2
#                 }
#             })
#
#
# for ds in ds_names:
#     print(ds)
#     pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + ds + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')  # Get the dataframe
#     pu['fold_growth'] = np.exp(pu['fold_growth'])
#     temp_dict = {}  # where we will save all the regression data
#     if ds == 'Pooled_SM':
#         for dataset in pu.dataset.unique():
#             ds_cond = (pu['dataset'] == dataset)
#             for lin_id in pu[ds_cond].lineage_ID.unique():
#                 relevant = pu[ds_cond & (pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
#
#                 if relevant.generation.max() < 29:  # We only want long enough traces
#                     continue
#                 else:
#                     print('lineage: {}, lin_id: {}, max_gen: {}'.format(dataset, lin_id, relevant.generation.max()))
#
#                 get_regression_data()
#     else:
#         for lin_id in pu.lineage_ID.unique():
#             relevant = pu[(pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
#
#             if relevant.generation.max() < 29:  # We only want long enough traces
#                 continue
#             else:
#                 print('lin_id: {}, max_gen: {}'.format(lin_id, relevant.generation.max()))
#
#             get_regression_data()
#
#     # Append them to the big dataframes
#     regression_dataframe = regression_dataframe.append(pd.DataFrame.from_dict(temp_dict, "index"), ignore_index=True)
#
# print('Finished!')
#
# regression_dataframe.to_csv(os.path.dirname(os.path.abspath(__file__))+'/regression_dataframe_exp_fg.csv', index=False)

regression_dataframe = pd.read_csv(os.path.dirname(os.path.abspath(__file__))+'/regression_dataframe_exp_fg.csv')

regression_dataframe = regression_dataframe[regression_dataframe['ds'].isin(ds_names)]  # Because of Maryam_LongTraces

save_folder = os.path.dirname(os.path.abspath(__file__))+'/Figures'
create_folder(save_folder)

print(regression_dataframe)
seaborn_preamble()
fig, ax = plt.subplots(tight_layout=True, figsize=[13, 6.5])
for var1, var2 in same_cell_variable_combos:
    s1 = symbols['physical_units'][var1] if var1 != 'fold_growth' else r'$e^{\phi}$'
    s2 = symbols['physical_units'][var2] if var2 != 'fold_growth' else r'$e^{\phi}$'
    cond = (regression_dataframe['var1'] == var1) & (regression_dataframe['var2'] == var2)
    sns.pointplot(data=regression_dataframe[cond], x='n', y='slope', hue='ds', palette=sns.color_palette("tab10"), dodge=True, linewidth=.1, capsize=.1, capthick=.1)   # , linestyles=['-', '--', ':', '-.']
    plt.ylabel('regression slope of {} {}'.format(s1, s2))
    plt.xlabel('')
    # plt.xticks(np.arange(min(regression_dataframe[cond].n), max(regression_dataframe[cond].n) + 1, 15), np.arange(min(regression_dataframe[cond].n), max(regression_dataframe[cond].n) + 1, 15))
    plt.legend(title='')
    plt.tight_layout()
    # plt.savefig('{}/{} {}.png'.format(save_folder, var1, var2), dpi=300)
    plt.show()
    plt.close()

# sns.scatterplot(data=relevant, x='length_birth', y='fold_growth')
# plt.show()
# plt.close()
