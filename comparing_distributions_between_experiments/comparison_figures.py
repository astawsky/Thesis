#!/usr/bin/env bash

import pandas as pd
import numpy as np
from CustomFuncsAndVars.global_variables import symbols, units, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns

# pu_of_datasets = {}

###############################################################################################################

seaborn_preamble()

for count, data_origin in enumerate(dataset_names):
    if data_origin == 'Pooled_SM':
        continue

    save_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin
    pu = pd.read_csv(save_folder+'/physical_units.csv')
    # pu_of_datasets.update({data_origin: pu})

    sns.kdeplot(data=pu, x='generationtime', y='growth_rate', color=cmap[count], label=data_origin)
plt.xlabel(symbols['physical_units']['generationtime'] + ' ' + units['generationtime'])
plt.ylabel(symbols['physical_units']['growth_rate'] + ' ' + units['growth_rate'])
plt.ylim(top=4, bottom=0)
plt.xlim(right=2.2)
plt.legend()
# plt.savefig('growth rate and generationtime.png', dpi=300)
plt.show()
plt.close()

###############################################################################################################

seaborn_preamble()

for count, data_origin in enumerate(dataset_names):
    if data_origin == 'Pooled_SM':
        continue

    save_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin
    pu = pd.read_csv(save_folder+'/physical_units.csv')
    # pu_of_datasets.update({data_origin: pu})

    sns.kdeplot(data=pu, x='length_birth', y='growth_rate', color=cmap[count], label=data_origin)
plt.xlabel(symbols['physical_units']['length_birth'] + ' ' + units['length_birth'])
plt.ylabel(symbols['physical_units']['growth_rate'] + ' ' + units['growth_rate'])
plt.ylim(top=4, bottom=0)
plt.xlim(left=0, right=7.2)
plt.legend()
# plt.savefig('growth rate and length_birth.png', dpi=300)
plt.show()
plt.close()

seaborn_preamble()

###############################################################################################################

for count, data_origin in enumerate(dataset_names):
    if data_origin == 'Pooled_SM':
        continue
    
    save_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin
    pu = pd.read_csv(save_folder+'/physical_units.csv')
    # pu_of_datasets.update({data_origin: pu})
    
    sns.kdeplot(data=pu, x='generationtime', y='length_birth', color=cmap[count], label=data_origin)
plt.xlabel(symbols['physical_units']['generationtime'] + ' ' + units['generationtime'])
plt.ylabel(symbols['physical_units']['length_birth'] + ' ' + units['length_birth'])
plt.ylim(bottom=0, top=7)
plt.xlim(left=0, right=4)
plt.legend()
plt.savefig('length_birth and generationtime.png', dpi=300)
# plt.show()
plt.close()

