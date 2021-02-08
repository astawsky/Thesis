#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names, shuffle_lineage_generations
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations
import os
from sklearn.neighbors import KernelDensity


def main():
    for ds in ['Pooled_SM']:
        print(ds)
        
        pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
        pu['experiment'] = ds
        pu['division_ratio'] = np.log(pu['division_ratio'])
        
        j = {len(pu[pu['lineage_ID'] == lin_id]): lin_id for lin_id in pu.lineage_ID.unique()}

        dim = 3
        projection = None if dim == 2 else '3d'
        
        for lin_id in pu.lineage_ID.unique():
            print(lin_id)
            
            print(np.mean(pu['division_ratio'] + pu['fold_growth']), np.std(pu['division_ratio'] + pu['fold_growth']))
            
            rel = pu[pu['lineage_ID'] == lin_id].copy().sort_values('generation')[['generationtime', 'growth_rate', 'division_ratio']].values

            fig = plt.figure(figsize=(12, 6))

            ax = fig.add_subplot(111, projection=projection)
            
            # print(*np.hsplit(rel, dim))
            # print(rel)

            ax.scatter(*np.hsplit(rel, dim), c=np.arange(rel.shape[0]), cmap=cm.binary)
            ax.plot3D(np.hsplit(rel, dim)[0].flatten(), np.hsplit(rel, dim)[1].flatten(), np.hsplit(rel, dim)[2].flatten())  #
            ax.set_xlabel(symbols['physical_units']['generationtime'])
            ax.set_ylabel(symbols['physical_units']['growth_rate'])
            ax.set_zlabel(r'log($f$)')  #symbols['physical_units']['division_ratio'])
            plt.show()
            plt.close()
        
        
        
        exit()
        
        print(j[np.max(list(j.keys()))], np.max(list(j.keys())))
        
        rel = pu[pu['lineage_ID'] == j[np.max(list(j.keys()))]].copy().sort_values('generation')[['generationtime', 'growth_rate', 'division_ratio']].values
        print(rel)

        
        


if __name__ == '__main__':
    main()
