#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations
import os

total_df = pd.DataFrame()
for data_origin in sm_datasets[:-1]:
    print(data_origin)
    pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/physical_units.csv'
    # pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'

    pu = pd.read_csv(pu)
    pu['experiment'] = data_origin

    total_df = total_df.append(pu, ignore_index=True)
    
pop_std = total_df[phenotypic_variables].std()

cv_df = pd.DataFrame(columns=['Ensemble', 'to_average', 'variable'])
for k, v in (total_df[phenotypic_variables].std()).to_dict().items():
    cv_df = cv_df.append({
        'Ensemble': 'Pooled',
        'to_average': 1,
        'variable': k
    }, ignore_index=True)

every_exp = total_df.experiment.unique()
for experiment in every_exp:
    print(experiment)

    exp_cond = (total_df['experiment'] == experiment)

    df = total_df[exp_cond].copy()

    for k, v in (df[phenotypic_variables].std()).to_dict().items():
        cv_df = cv_df.append({
            'Ensemble': 'Exp',
            'to_average': (len(every_exp) * v) / pop_std[k],
            'variable': k
        }, ignore_index=True)

    every_trap = df.trap_ID.unique()
    for trap_id in every_trap:
        trap_cond = (df['trap_ID'] == trap_id)

        df1 = df[trap_cond].copy()

        for k, v in (df1[phenotypic_variables].var()).to_dict().items():
            cv_df = cv_df.append({
                'Ensemble': 'Env',
                'to_average':  (len(every_trap) * v) / pop_std[k],
                'variable': k
            }, ignore_index=True)

        for trace in ['A', 'B']:
            trace_cond = (df1['trace'] == trace)

            df2 = df1[trace_cond].copy()

            for k, v in (df2[phenotypic_variables].std()).to_dict().items():
                cv_df = cv_df.append({
                    'Ensemble': 'Lin',
                    'to_average':  (len(every_trap) * v) / pop_std[k],
                    'variable': k
                }, ignore_index=True)

    # print(cv_df)

cv_df.to_csv('/Users/alestawsky/PycharmProjects/Thesis/cv_per_ensemble/cv_df.csv', index=False)

cv_df = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/cv_per_ensemble/cv_df.csv')

cv_df = cv_df.replace(symbols['physical_units'])

seaborn_preamble()
fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])

# sns.boxplot(data=cv_df, x='Ensemble', y='to_average', hue='variable', showfliers=False, order=['Lin', 'Env', 'Exp'], hue_order=[symbols['physical_units'][l] for l in phenotypic_variables])
sns.pointplot(data=cv_df, x='Ensemble', y='to_average', hue='variable', ci=None, dodge=True, order=['Lin', 'Env', 'Exp', 'Pooled'], hue_order=[symbols['physical_units'][l] for l in phenotypic_variables])
plt.xlabel('')
plt.legend(title='')
plt.savefig('pointplot.png', dpi=300)
plt.show()
plt.close()
