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

for ds in ['MG1655_inLB_LongTraces']:  # dataset_names:
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    
    for lin_id in pu.lineage_ID.unique():
        mask_lin_id = (pu['lineage_ID'] == lin_id)
        
        fig, axes = plt.subplots(2, 1)
        for var1, ax in zip(['generationtime', 'growth_rate', 'division_ratio', 'length_birth'], [axes[0], axes[0], axes[0], axes[1]]):  # fold_growth
            cumsum = np.cumsum((
                    (pu[mask_lin_id].sort_values('generation', ascending=0).dropna()[var1] - pu[mask_lin_id].sort_values('generation', ascending=0)[var1].dropna().mean()) /
                    pu[mask_lin_id].sort_values('generation', ascending=0)[var1].dropna().std()
            ).values)
            # ax.plot(pu[mask_lin_id].sort_values('generation').dropna()['generation'].values, cumsum, label=symbols['physical_units'][var1])
            ax.plot(pu[mask_lin_id].sort_values('generation').dropna()['generation'].values, cumsum, label=symbols['physical_units'][var1])
        
        plt.legend()
        plt.show()
        plt.close()
