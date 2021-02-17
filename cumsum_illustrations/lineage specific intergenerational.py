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
from matplotlib import cm

var1 = 'length_birth'
# var2 = var1
var2 = 'generationtime'

for ds in ['MG1655_inLB_LongTraces']:  # dataset_names:
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    
    fig, axes = plt.subplots(2, 1, tight_layout=True)
    axes[0].axhline(0, c='k')
    axes[1].axhline(0, c='k')
    
    for lin_id in np.random.choice(pu.lineage_ID.unique(), size=12, replace=False):
        rel = pu[pu['lineage_ID'] == lin_id].copy()[[var1, var2, 'generation']].dropna().sort_values('generation')
        rel = rel.loc[:, ~rel.columns.duplicated()]
        
        distances = np.arange(min(len(rel), 30))
        
        corr_per_dist = [pearsonr(rel[var1].values, rel[var2].values)[0] if distance == 0 else pearsonr(rel[var1].values[distance:], rel[var2].values[:-distance])[0] for distance in
                         distances]
        
        corr_per_dist_scrambled = [pearsonr(rel[var1].sample(frac=1).values, rel[var2].sample(frac=1).values)[0] if distance == 0 else
                                   pearsonr(rel[var1].sample(frac=1).values[distance:], rel[var2].sample(frac=1).values[:-distance])[0] for distance in
                                   distances]
        
        axes[0].plot(distances, corr_per_dist)
        axes[1].plot(distances, corr_per_dist_scrambled)
    axes[0].set_title('ordered')
    axes[1].set_title('shuffled')
    axes[0].set_ylim([-.4, 1])
    axes[1].set_ylim([-.4, 1])
    
    plt.show()
    plt.close()
