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

pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MG1655_inLB_LongTraces/ProcessedData/z_score_under_3/physical_units_without_outliers.csv').sort_values(['lineage_ID', 'generation'])
# pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv').sort_values(['lineage_ID', 'generation'])

scale = 1

# for lin_id in pu.lineage_ID.unique():
#
#     print(lin_id)

lin_id = 14
    
# stylistic reasons
sns.set_context('paper', font_scale=scale)
sns.set_style("ticks", {'axes.grid': True})

fig, axes = plt.subplots(2, 1, tight_layout=True, figsize=[6.5 * scale, 4.2 * scale])
axes[1].set_xlabel('generation')

for ax, variable in zip(axes, ['length_birth', 'growth_rate']):
    to_plot = pu[pu['lineage_ID'] == lin_id][variable].values
    
    print(len(to_plot))

    if len(to_plot) < 50:
        continue
        
    slope, intercept = linregress(np.arange(len(pu[pu['lineage_ID'] == lin_id][variable].dropna().values)),
                                  pu[pu['lineage_ID'] == lin_id][variable].dropna().values)[:2]
    ax.set_ylabel(symbols['physical_units'][variable] + ' ' + units[variable])
    ax.axhline(np.nanmean(to_plot), ls='-', c='k')
    ax.plot(to_plot, marker='o', color=cmap[0])
    ax.plot(intercept + (slope * np.arange(len(to_plot))), ls='--', color=cmap[1])
# plt.show()
# plt.savefig('Figures/dfa, trend illustration.png', dpi=300)
# plt.close()

