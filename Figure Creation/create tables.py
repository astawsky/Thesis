#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os

# growth_conditions_latex = pd.DataFrame(columns=['Exp.',
#                                                 'Date (d/m/y)',
#                                                 'Temp. (C)',
#                                                 'Medium',
#                                                 'Reference'])
#
# dates = ['15/10/18', '27/6/18', '13/7/18', '28/7/18', '12/10/18', 'Lambda LB ?', 'Maryam_LT (21/8/19, 25/9/19, ???)', 'MG1655 LT ?', '', '', '', '29/5/9', '30/9/9', '23/9/9', '22/9/9', '10/2/9',
#          '29/1/9', '2/7/9', '31/1/9', '25/5/9', '12/5/9']
#
# temp = ['', '', '', '', '', '?', '?', '?', '25', '27', '37', ]


references = [
    'Vashistha et al. 2021', 'Vashistha et al. 2021', 'Vashistha et al. 2021', 'Vashistha et al. 2021', 'Vashistha et al. 2021',
    'Susman et al. 2018', 'Korham et al. 2020', 'Susman et al. 2018', 'Tanouchi et al. 2015 (25C)', 'Tanouchi et al. 2015 (27C)', 'Tanouchi et al. 2015 (37C)',
    'Wang et al. 2010 (29/5/9)', 'Wang et al. 2010 (30/9/9)', 'Wang et al. 2010 (23/9/9)', 'Wang et al. 2010 (22/9/9)', 'Wang et al. 2010 (10/2/9)', 'Wang et al. 2010 (29/1/9)',
    'Wang et al. 2010 (2/7/9)', 'Wang et al. 2010 (31/1/9)', 'Wang et al. 2010 (25/5/9)', 'Wang et al. 2010 (12/5/9)'
]

gc = {count: {'Exp.': count, 'Ref.': ref} for ds, ref, count in zip(dataset_names[:-1], references, range(len(references))) if ds != 'Pooled_SM'}

growth_conditions = pd.DataFrame.from_dict(gc, "index")

print(growth_conditions)
print(dataset_names[6:-1])
print(growth_conditions.to_latex(float_format="{:0.2f}".format, index=False))

exit()

length_latex = pd.DataFrame(columns=['Exp.',
                                     '# of lineages',
                                     '# of cycles per lineage',
                                     '# lineages with at least 30 cycles'])
cv_latex = pd.DataFrame(columns=['Experiment Label'] + [symbols['physical_units'][p] for p in phenotypic_variables])

count = 0
for ds in dataset_names[:-1]:
    if ds == 'Pooled_SM':
        continue
    print(count, ds)
    
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    
    to_add = {'Experiment Label': count}
    
    to_add.update({symbols['physical_units'][p]: pu[p].std() / pu[p].mean() for p in phenotypic_variables})
    
    cv_latex = cv_latex.append(to_add, ignore_index=True)
    
    count += 1

# print(cv_latex.mean())
#
# print(cv_latex.std())
#
# print(cv_latex.to_latex(float_format="{:0.2f}".format, index=False))
