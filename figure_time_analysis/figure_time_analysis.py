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

pu = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/physical_units.csv')

sm = True

if sm:
    for dataset in pu.dataset.unique():
        ds_cond = (pu['dataset'] == dataset)
        for lin_id in pu[ds_cond].lineage_ID.unique():
            relevant = pu[ds_cond & (pu['lineage_ID'] == lin_id)].copy().sort_values('generation')
            
            if relevant.generation.max() < 20:  # We only want long enough traces
                continue
            print(relevant)
            
            sns.scatterplot(data=relevant, x='length_birth', y='fold_growth')
            plt.show()
            plt.close()
            
            exit()
else:
    pass
