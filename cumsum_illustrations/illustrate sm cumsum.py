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

for ds in ['Pooled_SM']:  # dataset_names:
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['division_ratio'] = np.log(pu['division_ratio'])
    pu = pu[pu['dataset'] == 'SL'].copy()

    for trap_id in pu.trap_ID.unique():
        print(trap_id)
        a_array = pu[(pu['trap_ID'] == trap_id) & (pu['trace'] == 'A')].copy().sort_values('generation')
        b_array = pu[(pu['trap_ID'] == trap_id) & (pu['trace'] == 'A')].copy().sort_values('generation')
        ax = plt.axes(projection='3d')
        ax.scatter(a_array.generationtime.values, a_array.growth_rate.values, a_array.division_ratio.values, c=np.arange(len(a_array)), cmap=cm.binary)
        ax.plot3D(a_array.generationtime.values, a_array.growth_rate.values, a_array.division_ratio.values, color='blue', label='A')
        ax.scatter(b_array.generationtime.values, b_array.growth_rate.values, b_array.division_ratio.values, c=np.arange(len(a_array)), cmap=cm.binary)
        ax.plot3D(b_array.generationtime.values, b_array.growth_rate.values, b_array.division_ratio.values, color='red', label='B')
        ax.set_xlabel(r'$\tau$')
        ax.set_ylabel(r'$\alpha')
        ax.set_zlabel(r'$f')
        plt.legend()
        plt.show()
        plt.close()

exit()

for ds in ['MG1655_inLB_LongTraces']:  # dataset_names:
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['division_ratio'] = np.log(pu['division_ratio'])
    
    for lin_id in pu.lineage_ID.unique():
        mask_lin_id = (pu['lineage_ID'] == lin_id)
        
        possible_gens = pu[mask_lin_id].sort_values('generation')['generation'].values
        fold_growths = pu[mask_lin_id].sort_values('generation')['fold_growth'].values
        division_ratios = pu[mask_lin_id].sort_values('generation')['division_ratio'].values
        grs = pu[mask_lin_id].sort_values('generation')['growth_rate'].values
        gts = pu[mask_lin_id].sort_values('generation')['generationtime'].values
        lbs = pu[mask_lin_id].sort_values('generation')['length_birth'].values

        fold_then_div = np.array([divs*np.exp(fold) if gen+1 in possible_gens else np.nan for gen, fold, divs in zip(possible_gens, fold_growths[:-1], division_ratios[1:])])
        grs = np.array([gr if gen+1 in possible_gens else np.nan for gen, gr in zip(possible_gens[1:], grs)])
        gts = np.array([gt if gen+1 in possible_gens else np.nan for gen, gt in zip(possible_gens[1:], gts)])
        
        normalize = lambda x: (x - np.nanmean(x)) / np.nanstd(x)
        
        # fig, axes = plt.subplots(2, 1)
        # axes[0].axhline(0, ls='-', c='k')
        # axes[1].axhline(0, ls='-', c='k')
        # repos = {}
        # for var1, ax in zip(['generationtime', 'growth_rate', 'division_ratio', 'length_birth', 'fold_growth', 'div_and_fold'], [axes[0], axes[0], axes[1], axes[0], axes[1], axes[1]]):  # fold_growth
        #     cumsum = np.cumsum((
        #                                (pu[mask_lin_id].sort_values('generation', ascending=1).dropna()[var1] - pu[mask_lin_id].sort_values('generation', ascending=1)[var1].dropna().mean()) /
        #                                pu[mask_lin_id].sort_values('generation', ascending=1)[var1].dropna().std()
        #                        ).values)
        #     # repos.update({var1: cumsum})
        #     # ax.plot(pu[mask_lin_id].sort_values('generation').dropna()['generation'].values, cumsum, label=symbols['physical_units'][var1])
        #     ax.plot(pu[mask_lin_id].sort_values('generation').dropna()['generation'].values, cumsum, label=symbols['physical_units'][var1])
        #
        # axes[0].legend()
        # axes[1].legend()
        # # plt.legend()
        # plt.show()
        # plt.close()
        
        fig, ax = plt.subplots()
        
        plt.plot(normalize(fold_then_div), color='blue')
        plt.plot(normalize(grs), color='orange')
        plt.plot(normalize(gts), color='green')
        plt.show()
        plt.close()
