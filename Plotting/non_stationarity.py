#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, mm_data_names, symbols, seaborn_preamble, shuffle_info, create_folder
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress


def ergodicity_breaking_in_lineage(lineage, min_ws, max_ws=None, window_size_steps=5):
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Maximum window size possible
    if not max_ws:
        max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we might skip some
    window_sizes = np.arange(min_ws, max_ws, window_size_steps, dtype=int)
    
    # This is where we keep the rescaled ranges for all partial time-series
    mse_array = []
    
    # if please:
    #     plt.plot(lineage)
    #     plt.show()
    #     plt.close()
    #
    #     plt.plot(total_walk)
    #     plt.show()
    #     plt.close()
    
    for ws in window_sizes:
        # The windowed/partial time-series
        partial_time_series = np.array([lineage.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, ws, dtype=int)[:-1]])
        
        print([[starting_point, starting_point + ws] for starting_point in np.arange(0, total_length, ws, dtype=int)[:-1]])
        print(lineage)
        print(partial_time_series)
        exit()
        
        # Where we will keep the mean squared error for each window we have
        window_fs = []
        
        # line time-series
        line_ts = np.array([np.array([np.mean(ts) for _ in np.arange(ts)]) for ts in partial_time_series])
        
        # Necessary
        assert len(line_ts) == len(partial_time_series)
        
        print()
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            
            # if np.isnan(ts).any():
            #     print('ts.isnull().values.any()')
            #     print(ts)
            #     # Get the linear regression parameters
            #     slope, intercept, _, _, _ = linregress(np.arange(len(ts)), ts)
            #     print(slope)
            #     print(intercept)
            #     exit()
            
            # Get the linear regression parameters
            slope, intercept, _, _, _ = linregress(np.arange(len(ts)), ts)
            
            # assert slope
            
            # Linear approximation
            line = [intercept + slope * dom for dom in np.arange(len(ts))]
            
            # Necessary
            assert len(ts) == len(line)
            
            # Mean Squared Deviation of the linear approximation
            f = np.sqrt(((ts - line) ** 2).sum() / len(ts))
            
            # Add it to the Mean Squared Error of the window size ws
            window_fs.append(f)
        
        # Add it to the Mean Squared Error of the window size ws
        mse_array.append(np.mean(window_fs))
    
    # So that we can regress successfully
    assert len(window_sizes) == len(mse_array)
    
    # if np.isnan(np.log(window_sizes)).any():
    #     print(window_sizes)
    #     print(np.log(window_sizes))
    #     exit()
    #
    # if np.isnan(np.log(mse_array)).any():
    #     print(mse_array)
    #     print(np.log(mse_array))
    #     exit()
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(window_sizes), np.log(mse_array))
    
    # # See what it looks like
    # plt.scatter(window_sizes, mse_array)
    # plt.plot(window_sizes, [np.exp(intercept) * (l ** slope) for l in window_sizes], label='{:.2}'.format(slope))
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.legend()
    # plt.show()
    # plt.close()
    # # exit()
    
    return [window_sizes, mse_array, slope, intercept, std_err]


def main(args):
    # import the labeled measured bacteria in physical units
    physical_units = pd.read_csv('{}/physical_units_with_outliers.csv'.format(args.save_folder))
    artificial_lineages = shuffle_info(physical_units, mm=True)
    
    seaborn_preamble()
    
    # This is to see the differences between the slope or intercepts of different lineages from different datasets and using different scaling analyses
    histogram_of_regression = pd.DataFrame(columns=['dataset', 'lineage_ID', 'slope', 'intercept', 'std_err', 'variable', 'kind'])
    
    # This is to recreate the loglog plots whose slopes give us the scalings of lineages based on some scaling analysis
    loglog_plot_recreation = pd.DataFrame(columns=['dataset', 'lineage_ID', 'window_sizes', 'y_to_regress', 'best_fit_line', 'variable', 'kind'])
    
    for variable in phenotypic_variables:
        
        print(variable)
        
        # fig, ax = plt.subplots(tight_layout=True)
        
        for lin_id in physical_units.lineage_ID.unique():
            print(lin_id)
            
            # The lineages that come from different "datasets"
            trace = physical_units[physical_units['lineage_ID'] == lin_id][variable].dropna()
            artificial = artificial_lineages[artificial_lineages['lineage_ID'] == lin_id][variable].dropna()  # .reset_indices(drop=True)
            white_noise = pd.Series(np.random.normal(0, 1, len(trace)))
            # noisy_wave = pd.Series([np.sin()])
            
            # What type of lineages do we want to look at?
            types_of_lineages, names_types_of_lineages = [trace, artificial, white_noise], ["Trace", "Artificial", "White Noise"]
            
            for lineage, dataset in zip(types_of_lineages, names_types_of_lineages):
    
                ergodicity_breaking_in_lineage(lineage, min_ws=3, max_ws=None, window_size_steps=2)
                exit()
                    
                # # Calculate the scaling of this "lineage" using this "kind" of analysis
                # window_sizes, y_to_regress, slope, intercept, std_err = scaling_analysis(lineage, min_ws=5, max_ws=None, window_size_steps=4, steps_between_windows=3)
                
                # # Save it
                # histogram_of_regression = histogram_of_regression.append(
                #     {'dataset': dataset, 'lineage_ID': lin_id, 'slope': slope, 'intercept': intercept, 'std_err': std_err, 'variable': symbols['physical_units'][variable], 'kind': kind},
                #     ignore_index=True)
                #
                # loglog_plot_recreation = loglog_plot_recreation.append(
                #     {
                #         'dataset': dataset, 'lineage_ID': lin_id, 'window_sizes': window_sizes, 'y_to_regress': y_to_regress,
                #         'best_fit_line': np.array([np.exp(intercept) * (ws ** slope) for ws in window_sizes]),
                #         'variable': symbols['physical_units'][variable], 'kind': kind
                #     }, ignore_index=True)
                #
                # # Make sure there are no NaNs
                # if histogram_of_regression.isnull().values.any():
                #     print(kind)
                #     print({'dataset': dataset, 'lineage_ID': lin_id, 'slope': slope, 'intercept': intercept, 'std_err': std_err, 'variable': symbols['physical_units'][variable], 'kind': kind})
                #     raise IOError('histogram_of_regression.isnull().values.any()')
                # if loglog_plot_recreation.isnull().values.any():
                #     print(kind)
                #     print({
                #         'dataset': dataset, 'lineage_ID': lin_id, 'window_sizes': window_sizes, 'y_to_regress': y_to_regress,
                #         'best_fit_line': np.array([np.exp(intercept) * (ws ** slope) for ws in window_sizes]),
                #         'variable': symbols['physical_units'][variable], 'kind': kind
                #     })
                #     raise IOError('loglog_plot_recreation.isnull().values.any()')
        
        print('*' * 200)
    
    print('scaling exponents', histogram_of_regression, sep='\n')
    print('recreate df', loglog_plot_recreation, sep='\n')
    
    # Save them
    histogram_of_regression.to_csv('{}/scaling_exponents.csv'.format(args.save_folder), index=False)
    loglog_plot_recreation.to_csv('{}/loglog_scaling_recreation.csv'.format(args.save_folder), index=False)
    
    print('finished saving them')


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in mm_data_names:
        print(data_origin)
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        main(args)
        
        # df = pd.read_csv('{}/scaling_exponents.csv'.format(args.save_folder))
        #
        # print(df)
        #
        # plot_histograms_of_scaling_exponents(args.figs_location, df)
        
        print('*' * 200)
        exit()
    
    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time
    
    data_origin = 'SM'
    
    print(data_origin)
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    args = parser.parse_args()
    
    main(args)
    
    # plot_histograms_of_scaling_exponents(args.figs_location, pd.read_csv('{}/scaling_exponents.csv'.format(args.save_folder)))
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
