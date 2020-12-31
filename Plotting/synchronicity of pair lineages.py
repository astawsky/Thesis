#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, dataset_names, symbols, seaborn_preamble, shuffle_info, create_folder
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress


def main():
    # import the labeled measured bacteria in physical units
    pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/SM/physical_units_with_control.csv')
    # pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/SM/trace_centered_with_control.csv')
    
    # Make sure the dataframe has no NaNs
    assert not pu['generationtime'].isnull().values.any()
    
    for dataset in pu.dataset.unique():
        print(dataset)
        
        tmp_df = pd.DataFrame(columns=['generation', 'difference'])
        for trap_id in pu[pu['dataset'] == dataset].trap_ID.unique():
            print(trap_id)
            
            lin_a = pu[(pu['dataset'] == dataset) & (pu['trap_ID'] == trap_id) & (pu['trace'] == 'A')].sort_values('generation')['generationtime'].copy().values
            lin_b = pu[(pu['dataset'] == dataset) & (pu['trap_ID'] == trap_id) & (pu['trace'] == 'B')].sort_values('generation')['generationtime'].copy().values
            
            # Make sure they are the same length
            assert len(lin_a) == len(lin_b)
            
            # add it to the dataframe
            tmp_df = tmp_df.append(pd.DataFrame({'generation': np.arange(1, len(lin_a)+1, dtype=int), 'difference': np.cumsum(lin_a) - np.cumsum(lin_b)}), ignore_index=True)
            
            print(tmp_df)
            exit()
            
        # for each generation get the variance of all the differences across pairs (130 for SL, 85 for NL and CTRL)
        


if __name__ == '__main__':
    main()
