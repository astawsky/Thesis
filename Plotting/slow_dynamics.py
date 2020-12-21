#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, mm_data_names, symbols, seaborn_preamble, shuffle_info
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress


def hurst(lineage):
    
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)

    # Minimum window size possible
    min_ws = 10
    max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, 1, dtype=int)
    # window_sizes = [int(np.floor(total_length / (2 ** l))) for l in np.arange(np.floor(np.log(total_length) / np.log(2)))]
    
    # print(window_sizes)
    # exit()
    
    # The regression of the log of this and the log of the window sizes is the Hurst exponent
    y_to_regress = []
    
    for ws in window_sizes:
        
        # The windowed/partial time-series
        partial_time_series = [lineage.iloc[starting_point:starting_point+ws].values for starting_point in np.arange(0, total_length, 2, dtype=int)[:-1]]
        # partial_time_series = [lineage.iloc[starting_point:starting_point+ws].values for starting_point in np.arange(0, total_length, ws, dtype=int)[:-1]]

        # This is where we keep the rescaled ranges for all partial time-series
        rescaled_ranges_array = []
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            # centered time-series
            y = ts - (ts.sum() / ws)
            
            # Their walk around the mean
            walk = y.cumsum()
            
            # The absolute range
            ts_range = walk.max() - walk.min()
            
            # For scaling
            ts_std = walk.std()
            
            if ts_std == 0:
                print('standard div of this partial time-series was zero, so we exclude it from the average')
                continue
            
            # Scaled absolute range
            ts_rescaled_range = ts_range / ts_std
            
            # add it to the array
            rescaled_ranges_array.append(ts_rescaled_range)
            
            # print(y, walk, ts_range, ts_std, ts_rescaled_range, sep='\n')
            # exit()
            
        # # Just a check... Unless the ts_std = 0
        # assert len(rescaled_ranges_array) == len(partial_time_series)
        
        # Get the average rescaled range
        average_rescaled_range = np.array(rescaled_ranges_array).sum() / len(rescaled_ranges_array)
        
        # Append it to the array
        y_to_regress.append(average_rescaled_range)
        
    # Convert it to a numpy array
    y_to_regress = np.array(y_to_regress)
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(window_sizes), np.log(y_to_regress))
    
    return [window_sizes, y_to_regress, slope, intercept, std_err]


def main(args):
    # import the labeled measured bacteria in physical units
    physical_units = pd.read_csv('{}/physical_units_with_outliers.csv'.format(args.save_folder))
    artificial_lineages = shuffle_info(physical_units, mm=True)

    seaborn_preamble()
    
    for variable in phenotypic_variables:
    
        print(variable)
        
        # Where we will keep the hurst exponents
        hurst_hist = {'Trace': [], 'Artificial': []}
        int_hist = {'Trace': [], 'Artificial': []}
        
        fig, ax = plt.subplots(tight_layout=True)
        
        for lin_id in physical_units.lineage_ID.unique():
            
            trace = physical_units[physical_units['lineage_ID'] == lin_id][variable].dropna()
            artificial = artificial_lineages[artificial_lineages['lineage_ID'] == lin_id][variable].dropna()#.reset_indices(drop=True)
            # artificial = pd.Series(np.random.normal(0, 1, len(trace)))
            
            window_sizes, y_to_regress, slope, intercept, std_err = hurst(trace)
            
            hurst_hist['Trace'].append(slope)
            int_hist['Trace'].append(intercept)

            plt.plot(window_sizes, y_to_regress, color='blue')
            # plt.plot(window_sizes, [np.exp(intercept) * (ws ** slope) for ws in window_sizes], ls='--', label=r'Trace: $H = {:.2} \pm {:.2}$'.format(slope, std_err))
            
            window_sizes, y_to_regress, slope, intercept, std_err = hurst(artificial)
            
            hurst_hist['Artificial'].append(slope)
            int_hist['Artificial'].append(intercept)
        
            plt.plot(window_sizes, y_to_regress, color='orange')
            # plt.plot(window_sizes, [np.exp(intercept) * (ws ** slope) for ws in window_sizes], ls='--', label=r'Art.: $H = {:.2} \pm {:.2}$'.format(slope, std_err))

        trace_avg_h = np.array(hurst_hist['Trace']).mean()
        artificial_avg_h = np.array(hurst_hist['Artificial']).mean()
        trace_avg_int = np.array(int_hist['Trace']).mean()
        artificial_avg_int = np.array(int_hist['Artificial']).mean()
        plt.plot(np.arange(3, 201), [np.exp(trace_avg_int) * (ws ** trace_avg_h) for ws in np.arange(3, 201)], ls='--', color='grey', lw=2)
        plt.plot(np.arange(3, 201), [np.exp(artificial_avg_int) * (ws ** artificial_avg_h) for ws in np.arange(3, 201)], ls='--', color='black', lw=2)
        plt.title(symbols['physical_units'][variable])
        plt.yscale('log')
        plt.xscale('log')
        ax.annotate('Trace: {:.2}\nArtificial: {:.2}'.format(trace_avg_h, artificial_avg_h), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
        # plt.legend()
        plt.show()
        plt.close()
            
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
