#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, shuffle_info, mm_data_names
import pandas as pd
import numpy as np
from hurst import compute_Hc, random_walk
import matplotlib.pyplot as plt


def main(args):
    # import the labeled measured bacteria in physical units
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    for variable in phenotypic_variables:
    
        print(variable)
        for lin_id in physical_units.lineage_ID.unique():
            
            lineage = physical_units[physical_units['lineage_ID'] == lin_id]
            
            if len(lineage) > 100:
                # Evaluate Hurst equation
                H, c, data = compute_Hc(lineage[variable], kind='price', simplified=True)
                
                # # Plot
                # f, ax = plt.subplots()
                # ax.plot(data[0], c * data[0] ** H, color="deepskyblue")
                # ax.scatter(data[0], data[1], color="purple")
                # ax.set_xscale('log')
                # ax.set_yscale('log')
                # ax.set_xlabel('Time interval')
                # ax.set_ylabel('R/S ratio')
                # ax.grid(True)
                # plt.show()
                
                print("H={:.4f}, c={:.4f}".format(H, c))
            
        print('*'*200)


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in mm_data_names:
        
        data_origin = 'MG1655_inLB_LongTraces'
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        # parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
        #                     required=False, default='population_lineages.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        # parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
        #                     required=False, default='ergodicity_breaking_parameter.csv')
        # parser.add_argument('-kld', '--kld', metavar='', type=str,
        #                     help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
        #                     required=False, default='kullback_leibler_divergences.csv')
        # parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        args = parser.parse_args()
        
        main(args)
        
        exit()
        print('*' * 200)
    
    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time
    
    data_origin = 'SM'
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                        required=False, default='population_lineages.csv')
    parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                        required=False, default='time_averages.csv')
    parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                        required=False, default='ergodicity_breaking_parameter.csv')
    parser.add_argument('-kld', '--kld', metavar='', type=str,
                        help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                        required=False, default='kullback_leibler_divergences.csv')
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    args = parser.parse_args()
    
    main(args)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
