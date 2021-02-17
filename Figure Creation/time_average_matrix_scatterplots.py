#!/usr/bin/env bash

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from AnalysisCode.global_variables import phenotypic_variables, symbols, seaborn_preamble, dataset_names, get_time_averages_df, sm_datasets, mm_datasets, pool_experiments


scale = 1.5

sns.set_context('paper', font_scale=scale/2)
sns.set_style("ticks", {'axes.grid': True})

fig, axes = plt.subplots(len(phenotypic_variables), len(phenotypic_variables), tight_layout=True, figsize=[6.5 * scale, 6.5 * scale])

pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
ta = get_time_averages_df(pu, phenotypic_variables).sort_values('lineage_ID')
ta = ta[ta['max_gen'] > 7][phenotypic_variables].drop_duplicates().reset_index(drop=True)
print(ta)

# mem = []
for row, var1 in enumerate(phenotypic_variables):
    print(var1)
    for col, var2 in enumerate(phenotypic_variables):
        # row = row
        ax = axes[row, col]
        if col >= row:
            ax.set_frame_on(False)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            if col != 0:
                axes[row, col].set_yticklabels('')
            if row != len(phenotypic_variables)-1:
                axes[row, col].set_xticklabels('')
                
            # sns.regplot(line_kws=)
            intercept, slope = linregress(ta[var1].values, ta[var2].values)[:2]
            ax.scatter(ta[var1].values, ta[var2].values)#, marker='o')
            ax.plot(np.linspace(np.nanmin(ta[var1]), np.nanmax(ta[var1])), intercept + np.linspace(np.nanmin(ta[var1]), np.nanmax(ta[var1])) * slope, ls='--', c='k')
            
    # mem.append(var1)
    
for ax, variable in zip(axes[len(phenotypic_variables)-1, :], phenotypic_variables):
    ax.set_xlabel(symbols['time_averages'][variable], size='xx-large')

for ax, variable in zip(axes[:, 0], phenotypic_variables):
    ax.set_ylabel(symbols['time_averages'][variable], size='xx-large')
    
plt.show()
plt.close()
