#!/usr/bin/env bash

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from AnalysisCode.global_variables import phenotypic_variables, symbols, seaborn_preamble, dataset_names, get_time_averages_df, sm_datasets, mm_datasets, pool_experiments


def pyramid_of_pairwise_covariances(pu, ta, tc, fig, axes, variables=phenotypic_variables, figurename='covariance decomposition, main variables', annot=True):
    def normalize_correlation(df, variables):
        cov = df.cov()
        corr_df = pd.DataFrame(columns=variables, index=variables, dtype=float)
        for param_1 in variables:
            for param_2 in variables:
                normalization = (pu[param_1].std() * pu[param_2].std())
                
                # # If the value lies in the noise range then make it zero
                # if -.1 < cov.loc[param_1, param_2] / normalization < .1:
                #     corr_df.loc[param_1, param_2] = float('nan')  # 0
                # else:
                #     corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
                corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
        
        return corr_df
    
    mask = np.ones_like(normalize_correlation(pu, variables))
    mask[np.tril_indices_from(mask)] = False
    vmax, vmin = 1, -1
    
    npu = normalize_correlation(pu, variables).rename(columns=symbols['physical_units'], index=symbols['physical_units'])
    nta = normalize_correlation(ta, variables).rename(columns=symbols['time_averages'], index=symbols['time_averages'])
    ntc = normalize_correlation(tc, variables).rename(columns=symbols['trace_centered'], index=symbols['trace_centered'])
    
    sns.heatmap(npu, annot=annot, center=0, vmax=vmax,
                vmin=vmin, cbar=False, ax=axes[0], mask=mask, square=True, fmt='.2f')
    axes[0].set_title('A', x=-.2, fontsize='xx-large')
    
    sns.heatmap(nta, annot=annot, center=0, vmax=vmax, vmin=vmin,
                cbar=False, ax=axes[1], mask=mask, square=True, fmt='.2f')
    axes[1].set_title('B', x=-.2, fontsize='xx-large')
    
    cbar_ax = fig.add_axes([.91, .1, .03, .8])
    
    sns.heatmap(ntc, annot=annot, center=0, vmax=vmax,
                vmin=vmin, cbar=True, ax=axes[2], mask=mask, cbar_kws={"orientation": "vertical"}, square=True, cbar_ax=cbar_ax, fmt='.2f')
    axes[2].set_title('C', x=-.2, fontsize='xx-large')
    
    # plt.suptitle(args.data_origin)
    
    fig.tight_layout(rect=[0, 0, .9, 1])
    
    plt.savefig('Figures/{}.png'.format(figurename), dpi=300)
    # plt.show()
    plt.close()

    
# The variables we want to plot
variables = ['div_and_fold', 'fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']

group = [sm_datasets[-1]]  # sm_datasets[:-1] + mm_datasets
name = 'Sister Machine Experiments'  # 'Hanna Experiments in LB'

# pu = pool_experiments(group=group, name=name, outliers=False, dimensions='pu')
# tc = pool_experiments(group=group, name=name, outliers=False, dimensions='tc')
# ta = pool_experiments(group=group, name=name, outliers=False, dimensions='ta')

pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
ta = get_time_averages_df(pu, phenotypic_variables)
tc = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/trace_centered_without_outliers.csv')

scale = 1.5

sns.set_context('paper', font_scale=.7 * scale)
sns.set_style("ticks", {'axes.grid': False})

fig, axes = plt.subplots(1, 3, figsize=[6.5 * scale, 2.1 * scale])

pyramid_of_pairwise_covariances(pu, ta, tc, fig, axes, variables=variables, figurename='covariance decomposition', annot=True)

