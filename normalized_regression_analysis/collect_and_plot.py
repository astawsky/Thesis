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

ds_names = [sm_datasets[-1]] + cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets

same_cell_variable_combos = list(combinations(['length_birth', 'generationtime', 'growth_rate', 'division_ratio'], 2))
same_cell_variable_combos = same_cell_variable_combos + [('fold_growth', 'length_birth'), ('length_final', 'length_birth'), ('fold_growth', 'division_ratio'), ('added_length', 'division_ratio')]

regression_dataframe = pd.DataFrame()

for ds in ds_names:
    print(ds)
    pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + ds + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')  # Get the dataframe
    temp_dict = {}
    if ds == 'Pooled_SM':
        for dataset in pu.dataset.unique():
            ds_cond = (pu['dataset'] == dataset)
            for lin_id in pu[ds_cond].lineage_ID.unique():
                relevant = pu[ds_cond & (pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
                
                if relevant.generation.max() < 29:  # We only want long enough traces
                    continue
                else:
                    print('ds: {}, lin_id: {}, max_gen: {}'.format(ds, lin_id, relevant.generation.max()))

                for gen in np.arange(29, relevant.generation.max(), 5):
                    print(gen)
                    gen_cond = (relevant['generation'] <= gen)
                    for var1, var2 in same_cell_variable_combos:
                        
                        x, y = relevant[gen_cond][var1], relevant[gen_cond][var2]
                        x, y = ((x-x.mean())/x.std()).values, ((y-y.mean())/y.std()).values
                        
                        print(len(x), len(y))
                        exit()
                    
                        slope, intercept, r_value, _, _ = linregress(x, y)
                        
                        temp_dict.update({len(temp_dict): {
                            'slope': slope,
                            'intercept': intercept,
                            'r_value': r_value,
                            'lin_id': lin_id,
                            'ds': dataset,
                            'n': gen+1
                        }})
            
                sns.scatterplot(data=relevant, x='length_birth', y='fold_growth')
                plt.show()
                plt.close()
            
                exit()
    else:
        pass
