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

ds_names = cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets  # sm_datasets[:-1]

for ds in ['Pooled_SM']:#ds_names:
    print(ds)
    for var1, var2 in [('length_birth', 'length_final')]:  # [('length_birth', 'growth_rate'), ('length_birth', 'generationtime')]:  #('length_birth', 'fold_growth'), ('generationtime', 'growth_rate')
        pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + ds + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')  # Get the dataframe
        # pu[var1] = np.log(pu[var1])
        # pu[phenotypic_variables] = pu[
        #     pu[phenotypic_variables] - pu[phenotypic_variables].mean() < (3*pu[phenotypic_variables].std())
        # ][phenotypic_variables]
        x, y = pu[var1].values, pu[var2].values
        sns.kdeplot(x, y, color='gray')
        fake_x = np.linspace(np.nanmin(x), np.nanmax(x))
    
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
                    
                    # plt.plot(fake_x, [intercept + l * slope for l in fake_x])
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
                # print(slope, intercept)
    
                # plt.plot(fake_x, [intercept + l * slope for l in fake_x])
    
        # plt.scatter(pu.length_birth.mean(), pu.fold_growth.mean(), marker='o', s=50, c='k', zorder=500)
        
        # plt.title(ds)
        plt.scatter([pu[pu['lineage_ID'] == lin_id][var1].mean() for lin_id in pu.lineage_ID.unique()],
                    [pu[pu['lineage_ID'] == lin_id][var2].mean() for lin_id in pu.lineage_ID.unique()], marker='x', c='black', zorder=500)
        plt.xlabel(symbols['physical_units'][var1])
        plt.ylabel(symbols['physical_units'][var2])
        plt.tight_layout()
        # plt.savefig(f'normalized_regression_analysis/{var1} {var2}.png', dpi=300)
        plt.show()
        plt.close()
