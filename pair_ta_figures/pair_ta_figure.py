#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, cut_uneven_pairs, add_control, add_control_and_cut_extra_intervals
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, linregress
import os


""" Returns dataframe with same size containing the time-averages of each phenotypic variable instead of the local value """


def get_time_averages_sm(info, variables):
    
    # We keep the trap means here
    means_df = pd.DataFrame()
    
    # specify a lineage
    for lin_id in info['lineage_ID'].unique():
        # print(lin_id)
        
        # the values of the lineage we get from physical units
        lineage = info[(info['lineage_ID'] == lin_id)].copy()
        
        # add its time-average
        to_add = {
            'trace': lineage.trace.unique()[0],
            'lineage_ID': lin_id,
            'max_gen': len(lineage)
        }
        to_add.update({param: np.mean(lineage[param]) for param in variables})
        to_add = pd.DataFrame(to_add, index=[0])
        means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    
    return means_df.reset_index(drop=True)


def plot(df, label, var1, var2):
    
    variables = [var1, var2]
    
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    
    latex_symbols = [symbols[label][v] for v in variables]
    unit_symbols = [units[v] for v in variables]
    
    # # Take out the outliers
    # df[variables] = df[variables].where(
    #     np.abs(df[variables] - df[variables].mean()) < (3 * df[variables].std()),
    #     other=np.nan
    # )
    
    # stylistic reasons
    sns.set_context('paper', font_scale=1.2)
    sns.set_style("ticks", {'axes.grid': True})
    
    fig, axes = plt.subplots(figsize=[5, 5], tight_layout=True)
    
    sym1 = latex_symbols[0]+r'$_A$'
    sym2 = latex_symbols[1]+r'$_B$'
    
    relevant = pd.DataFrame()
    # print(df[(df['trace'] == 'A')].sort_values(['experiment', 'lineage_ID'])[var1].values)
    relevant[sym1] = df[(df['trace'] == 'A')].sort_values(['experiment', 'lineage_ID'])[var1].values
    relevant[sym2] = df[(df['trace'] == 'B')].sort_values(['experiment', 'lineage_ID'])[var2].values
    relevant['experiment'] = df[(df['trace'] == 'B')].sort_values(['experiment', 'lineage_ID'])['experiment'].values
    
    relevant = relevant.dropna()
    
    pcorr = str(pearsonr(relevant[sym2].values, relevant[sym1].values)[0])[:4]
    slope, intercept, r_value, _, std_err = linregress(relevant[sym2], relevant[sym1])
    
    sns.regplot(data=relevant, x=sym2, y=sym1, line_kws={'color': 'black', 'ls': '--', 'lw': 2})
    # sns.scatterplot(data=relevant, x=sym2, y=sym1, hue='experiment')
    # plt.legend(title='')
    
    # sns.regplot(x=df[sym2], y=df[sym1], data=df, ax=axes, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
    
    axes.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(str(slope)[:4]), xy=(.5, .92), xycoords=axes.transAxes, fontsize=13, ha='center', va='bottom',
                  color='red',
                  bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    axes.set_ylabel(sym1 + r' ' + unit_symbols[0])
    axes.set_xlabel(sym2 + r' ' + unit_symbols[1])
    
    # plt.savefig(r'/Users/alestawsky/PycharmProjects/Thesis/sm_ta_correlations_figure/experiment_centered_gr.png', dpi=300)
    plt.show()
    plt.close()


# The variables we want to plot
main_variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']

p = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/physical_units.csv')

total_pu = pd.DataFrame()
total_tc = pd.DataFrame()
total_ta = pd.DataFrame()
for data_origin in sm_datasets[:-1]:
    print(data_origin)
    pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/physical_units.csv')
    tc = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/trace_centered.csv')
    
    ta = get_time_averages_sm(pu, phenotypic_variables).drop_duplicates().sort_values(['lineage_ID'])
    
    pool = pu[phenotypic_variables].mean()
    
    pu['experiment'] = data_origin
    pu[phenotypic_variables] = pu[phenotypic_variables] - pool
    total_pu = total_pu.append(pu, ignore_index=True)
    tc['experiment'] = data_origin
    tc[phenotypic_variables] = tc[phenotypic_variables]
    total_tc = total_tc.append(tc, ignore_index=True)
    ta['experiment'] = data_origin
    ta[phenotypic_variables] = ta[phenotypic_variables] - pool
    total_ta = total_ta.append(ta, ignore_index=True)

plot(total_ta, 'time_averages', 'division_ratio', 'fold_growth')
