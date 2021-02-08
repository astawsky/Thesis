#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names, shuffle_lineage_generations
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations
import os


def dfa_short(lineage, min_ws, max_ws=None, window_size_steps=5, steps_between_windows=2):
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Maximum window size possible
    if not max_ws:
        max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, window_size_steps, dtype=int)
    
    # Convert it into a random walk of some persistence
    total_walk = (lineage - lineage.mean()).cumsum()
    
    # This is where we keep the rescaled ranges for all partial time-series
    mse_array = []
    
    for ws in window_sizes:
        # The windowed/partial time-series
        partial_time_series = [total_walk.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, steps_between_windows, dtype=int)[:-1]]
        
        # Where we will keep the mean squared error for each window we have
        window_fs = []
        
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
        
        plt.show()
        plt.close()
    
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
    # plt.tight_layout()
    # plt.legend()
    # plt.show()
    # plt.close()
    # # exit()
    
    return [window_sizes, mse_array, slope, intercept, std_err]


def mean_short(lineage, min_ws, max_ws=None, window_size_steps=5, steps_between_windows=2):
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Maximum window size possible
    if not max_ws:
        max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, window_size_steps, dtype=int)
    
    # Convert it into a random walk of some persistence
    total_walk = (lineage - lineage.mean()).cumsum()
    
    # This is where we keep the rescaled ranges for all partial time-series
    mse_array = []
    
    for ws in window_sizes:
        # The windowed/partial time-series
        partial_time_series = [total_walk.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, steps_between_windows, dtype=int)[:-1]]
        
        # Where we will keep the mean squared error for each window we have
        window_fs = []
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            
            # Get the linear regression parameters
            slope, intercept, _, _, _ = linregress(np.arange(len(ts)), ts)
            
            # 0-order approximation
            line = [np.mean(ts) for _ in np.arange(len(ts))]
            
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
    # exit()
    
    return [window_sizes, mse_array, slope, intercept, std_err]


def plot_histograms_of_scaling_exponents(figure_folder, histogram_of_regression):
    
    # Histogram of H of every kind
    for kind in histogram_of_regression.kind.unique():
        seaborn_preamble()
        fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
        plt.axhline(.5, ls='-', c='k', linewidth=1, zorder=1)
        
        sns.pointplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=[symbols['physical_units'][v] for v in phenotypic_variables],
                    join=False, dodge=True, capsize=.5)  # , capthick=1
        # sns.boxplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=[symbols['physical_units'][v] for v in phenotypic_variables],
        #             showfliers=False)
        plt.legend(title='')
        # plt.title(kind)
        plt.ylabel('')
        plt.xlabel('')
        plt.savefig('{}/{}.png'.format(figure_folder, kind), dpi=300)
        # plt.show()
        plt.close()


def recreate_loglog_plots(loglog_plot_recreation, figure_folder, individuals=False):
    # Create the folders where we will save the figures
    if individuals:
        create_folder('{}/Individual Scalings'.format(figure_folder))
        for kind in loglog_plot_recreation.kind.unique():
            create_folder('{}/Individual Scalings/{}'.format(figure_folder, kind))
            for dataset in loglog_plot_recreation.dataset.unique():
                create_folder('{}/Individual Scalings/{}/{}'.format(figure_folder, kind, dataset))
    else:
        for kind in loglog_plot_recreation.kind.unique():
            create_folder('{}/{}'.format(figure_folder, kind))
    
    # For every kind of loglog plot
    for kind in loglog_plot_recreation.kind.unique():
        print(kind)
        for variable in loglog_plot_recreation.variable.unique():
            if individuals:
                for dataset in loglog_plot_recreation.dataset.unique():
                    for lin_id in loglog_plot_recreation.lineage_ID.unique():
                        relevant = loglog_plot_recreation[
                            (loglog_plot_recreation['kind'] == kind) & (loglog_plot_recreation['lineage_ID'] == lin_id) & (loglog_plot_recreation['variable'] == variable) & (
                                    loglog_plot_recreation['dataset'] == dataset)].copy()
                        
                        fig, ax = plt.subplots(tight_layout=True)
                        plt.scatter(relevant['window_sizes'].values, relevant['y_to_regress'].values)
                        plt.plot(relevant['window_sizes'].values, relevant['best_fit_line'].values,
                                 label=r'log(y)$= {:.2} + {:.2} \,$log(x)'.format(relevant.intercept.unique()[0], relevant.slope.unique()[0]))
                        plt.ylabel('Variation')
                        plt.xlabel('Window Size')
                        plt.yscale('log')
                        plt.xscale('log')
                        plt.legend(title='')
                        plt.tight_layout()
                        plt.show()
                        # plt.savefig('{}/Individual Scalings/{}/{}/{}.png'.format(figure_folder, kind, dataset, lin_id), dpi=300)
                        plt.close()
                        exit()
            else:
                relevant = loglog_plot_recreation[
                    (loglog_plot_recreation['kind'] == kind) & (loglog_plot_recreation['variable'] == variable) & (loglog_plot_recreation['dataset'].isin(['Trace', 'Shuffled']))].copy()
                
                seaborn_preamble()
                
                set_trace_legend, set_shuffled_legend = False, False
                
                wss_array = {'Trace': np.array([]), 'Shuffled': np.array([])}
                ytr_array = wss_array.copy()
                
                # For each lineage separately
                for ind in relevant.index:
                    # Convert the window sizes and mean squared error per windows size from the dataframe from a string to an array of integers or floats
                    wss = relevant['window_sizes'].loc[ind].strip('][').split(' ')
                    wss = [int(r.split('\n')[0]) for r in wss if r != '']
                    ytr = relevant['y_to_regress'].loc[ind].strip('][').split(', ')
                    ytr = [float(r) for r in ytr if r != '']
                    
                    # If it is a trace lineage plot it and include it in the legend
                    if relevant.loc[ind]['dataset'] == 'Trace':
                        # The constant shuffled color
                        color = cmap[0]
                        
                        # Add the analysis curve to the array of its dataset for the population regression
                        wss_array['Trace'] = np.append(wss_array['Trace'], wss)
                        ytr_array['Trace'] = np.append(ytr_array['Trace'], ytr)
                        
                        # Include it in the legend but not more than once
                        if set_trace_legend:
                            plt.plot(wss, ytr, color=color, marker='x')
                        else:
                            plt.plot(wss, ytr, color=color, marker='x', label='Trace')
                            set_trace_legend = True
                    # If it is a shuffled lineage plot it and include it in the legend
                    elif relevant.loc[ind]['dataset'] == 'Shuffled':
                        # The constant shuffled color
                        color = cmap[1]
                        
                        # Add the analysis curve to the array of its dataset for the population regression
                        wss_array['Shuffled'] = np.append(wss_array['Shuffled'], wss)
                        ytr_array['Shuffled'] = np.append(ytr_array['Shuffled'], ytr)
                        
                        # Include it in the legend but not more than once
                        if set_shuffled_legend:
                            plt.plot(wss, ytr, color=color, marker='x')
                            pass
                        else:
                            plt.plot(wss, ytr, color=color, marker='x', label='Shuffled')
                            set_shuffled_legend = True
                    # We do not want to see the white noise
                    else:
                        continue
                
                # Get the linear regression of all the trajectories pooled for each dataset
                slope_trace, intercept_trace, _, _, std_err_trace = linregress(np.log(np.array(wss_array['Trace']).flatten()), np.log(np.array(ytr_array['Trace']).flatten()))
                slope_art, intercept_art, _, _, std_err_art = linregress(np.log(np.array(wss_array['Shuffled']).flatten()), np.log(np.array(ytr_array['Shuffled']).flatten()))
                
                # Plot the best fit line and its parameters
                plt.plot(np.unique(wss_array['Trace']), [np.exp(intercept_trace) * (l ** slope_trace) for l in np.unique(wss_array['Trace'])], ls='--', color='blue', linewidth=3,
                         label=r'$' + str(np.round(intercept_trace, 2)) + r'n^{' + str(np.round(slope_trace, 2)) + r'\pm' + str(np.round(std_err_trace, 3)) + r'}$')
                plt.plot(np.unique(wss_array['Shuffled']), [np.exp(intercept_art) * (l ** slope_art) for l in np.unique(wss_array['Shuffled'])], ls='--', color='red', linewidth=3,
                         label=r'$' + str(np.round(intercept_art, 2)) + r'n^{' + str(np.round(slope_art, 2)) + r'\pm' + str(np.round(std_err_art, 3)) + r'}$')
                
                plt.title(variable)
                plt.ylabel(r'$F(n)$')
                plt.xlabel('n')
                plt.yscale('log')
                plt.xscale('log')
                plt.legend(title='')
                plt.tight_layout()
                # plt.show()
                plt.savefig('{}/{}/{}.png'.format(figure_folder, kind, variable), dpi=300)
                plt.close()
    # plt.scatter(window_sizes, y_to_regress, color='blue')


def main(args):
    
    min_ws = 5
    max_ws = None
    window_size_steps = 2
    steps_between_windows = 3
    
    # import the labeled measured bacteria in physical units
    physical_units = pd.read_csv(args['pu'])
    shuffled_lineages = shuffle_lineage_generations(physical_units, mm=True if args['data_origin'] not in sm_datasets else False)
    
    # This is to see the differences between the slope or intercepts of different lineages from different datasets and using different scaling analyses
    histogram_of_regression = pd.DataFrame(columns=['dataset', 'lineage_ID', 'slope', 'intercept', 'std_err', 'variable', 'kind'])
    
    # This is to recreate the loglog plots whose slopes give us the scalings of lineages based on some scaling analysis
    loglog_plot_recreation = pd.DataFrame(columns=['dataset', 'lineage_ID', 'window_sizes', 'y_to_regress', 'best_fit_line', 'variable', 'kind'])
    
    for variable in phenotypic_variables:
        
        print(variable)
        
        # fig, ax = plt.subplots(tight_layout=True)
        
        for lin_id in physical_units.lineage_ID.unique():
            # print(lin_id)
            
            # The lineages that come from different "datasets"
            trace = physical_units[physical_units['lineage_ID'] == lin_id][variable].dropna()
            shuffled = shuffled_lineages[shuffled_lineages['lineage_ID'] == lin_id][variable].dropna()  # .reset_indices(drop=True)
            # white_noise = pd.Series(np.random.normal(0, 1, len(trace)))
            
            # If we have enough points to get a slope using this lineage
            if len(np.arange(min_ws, len(trace) / 2, window_size_steps, dtype=int)) < 3:
                # print('not enough cycles')
                continue
            
            # What type of lineages do we want to look at?
            types_of_lineages, names_types_of_lineages = [trace, shuffled], ["Trace", "Shuffled"]
            
            for lineage, dataset in zip(types_of_lineages, names_types_of_lineages):
                for scaling_analysis, kind in zip([dfa_short, mean_short], ['dfa (short)', 'mean (short)']):
                    
                    # Calculate the scaling of this "lineage" using this "kind" of analysis
                    window_sizes, y_to_regress, slope, intercept, std_err = scaling_analysis(lineage, min_ws=min_ws, max_ws=max_ws, window_size_steps=window_size_steps,
                                                                                             steps_between_windows=steps_between_windows)
                    
                    # Save it
                    histogram_of_regression = histogram_of_regression.append(
                        {'dataset': dataset, 'lineage_ID': lin_id, 'slope': slope, 'intercept': intercept, 'std_err': std_err, 'variable': symbols['physical_units'][variable], 'kind': kind},
                        ignore_index=True)
                    
                    loglog_plot_recreation = loglog_plot_recreation.append(
                        {
                            'dataset': dataset, 'lineage_ID': lin_id, 'window_sizes': window_sizes, 'y_to_regress': y_to_regress,
                            'best_fit_line': np.array([np.exp(intercept) * (ws ** slope) for ws in window_sizes]),
                            'variable': symbols['physical_units'][variable], 'kind': kind
                        }, ignore_index=True)
                    
                    # Make sure there are no NaNs
                    if histogram_of_regression.isnull().values.any():
                        print(kind)
                        print({'dataset': dataset, 'lineage_ID': lin_id, 'slope': slope, 'intercept': intercept, 'std_err': std_err, 'variable': symbols['physical_units'][variable], 'kind': kind})
                        raise IOError('histogram_of_regression.isnull().values.any()')
                    if loglog_plot_recreation.isnull().values.any():
                        print(kind)
                        print({
                            'dataset': dataset, 'lineage_ID': lin_id, 'window_sizes': window_sizes, 'y_to_regress': y_to_regress,
                            'best_fit_line': np.array([np.exp(intercept) * (ws ** slope) for ws in window_sizes]),
                            'variable': symbols['physical_units'][variable], 'kind': kind
                        })
                        raise IOError('loglog_plot_recreation.isnull().values.any()')
        # ax.annotate('Trace: {:.2}\nShuffled: {:.2}'.format(trace_avg_h, shuffled_avg_h), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
        #             bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
        
        print('*' * 200)
    
    print('scaling exponents', histogram_of_regression, sep='\n')
    print('recreate df', loglog_plot_recreation, sep='\n')
    
    # Save them
    histogram_of_regression.to_csv('{}/scaling_exponents.csv'.format(args['Scaling_Exponents']), index=False)
    loglog_plot_recreation.to_csv('{}/loglog_scaling_recreation.csv'.format(args['LogLog_Recreation']), index=False)
    
    print('finished saving them')
    
    # # Plot the figures
    # plot_histograms_of_scaling_exponents(args['save_figs'], histogram_of_regression)
    # recreate_loglog_plots(loglog_plot_recreation, args['save_figs'], individuals=True)


if __name__ == '__main__':
    import argparse
    
    # The variables we want to plot
    main_variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']
    
    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')
    
    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    parser.add_argument('-kinds_of_correlations', '--kinds_of_correlations', metavar='', nargs="+", help='Calculate pearson and/or variance decomposition correlation?', required=False,
                        default=['decomposition', 'pearson'])
    parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
                        default=dict(zip(['phenotypic_variables'], [phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    # Do all the Mother and Sister Machine data
    for data_origin in ['lambda_LB']:  # input_args.dataset_names:
        print(data_origin)
        
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # create_folder(current_dir + '/Dataframes')
        # create_folder(current_dir + '/Dataframes/' + data_origin)
        
        create_folder(current_dir + '/Scaling Exponents')
        create_folder(current_dir + '/LogLog Recreation')
        create_folder(current_dir + '/Scaling Exponents/'+data_origin)
        create_folder(current_dir + '/LogLog Recreation/'+data_origin)
        
        processed_data = os.path.dirname(current_dir) + '/Datasets/' + data_origin + '/ProcessedData/'
        
        """
        data_origin ==> Name of the dataset we are analysing
        raw_data ==> Where the folder containing the raw data for this dataset is
        processed_data ==> The folder we will put the processed data in
        """
        
        args = {
            'data_origin': data_origin,
            # Data singularities, long traces with significant filamentation, sudden drop-offs
            'Scaling_Exponents': current_dir + '/Scaling Exponents/'+data_origin,
            'LogLog_Recreation': current_dir + '/LogLog Recreation/'+data_origin,
            # 'dataframes': current_dir + '/Dataframes/' + data_origin,
            'pu': processed_data + '/physical_units.csv'
        }
        
        main(args)
        
        df = pd.read_csv('{}/scaling_exponents.csv'.format(args['Scaling_Exponents']))
        
        plot_histograms_of_scaling_exponents(args['Scaling_Exponents'], df)
        
        df = pd.read_csv('{}/loglog_scaling_recreation.csv'.format(args['LogLog_Recreation']))
        
        recreate_loglog_plots(df, args['LogLog_Recreation'], individuals=False)
