#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, mm_data_names, symbols, seaborn_preamble, shuffle_info, create_folder, cmap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from ast import literal_eval
from json import loads


def hurst_original(lineage, min_ws):
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Minimum window size possible
    max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, 1, dtype=int)
    # window_sizes = [int(np.floor(total_length / (2 ** l))) for l in np.arange(np.floor(np.log(total_length) / np.log(2)))]
    
    # print(window_sizes)
    
    # The regression of the log of this and the log of the window sizes is the Hurst exponent
    y_to_regress = []
    
    for ws in window_sizes:
        
        # The windowed/partial time-series
        partial_time_series = [lineage.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, 2, dtype=int)[:-1]]
        partial_time_series = [lineage.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, ws, dtype=int)[:-1]]
        
        # This is where we keep the rescaled ranges for all partial time-series
        rescaled_ranges_array = []
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            # centered time-series
            y = ts - (ts.sum() / len(ts))
            
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
    
    # print(window_sizes)
    # print(y_to_regress)
    # print(len(window_sizes))
    # print(len(y_to_regress))
    # exit()
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(window_sizes), np.log(y_to_regress))
    
    return [window_sizes, y_to_regress, slope, intercept, std_err]


def dfa_original(lineage, min_ws):
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Minimum window size possible
    max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, 1, dtype=int)
    # window_sizes = [int(np.floor(total_length / (2 ** l))) for l in np.arange(np.floor(np.log(total_length) / np.log(2)))]
    
    # print(window_sizes)
    
    # Convert it into a random walk of some persistence
    total_walk = (lineage - lineage.mean()).cumsum()
    
    # This is where we keep the rescaled ranges for all partial time-series
    mse_array = []
    
    for ws in window_sizes:
        # The windowed/partial time-series
        # partial_time_series = [lineage.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, 2, dtype=int)[:-1]]
        partial_time_series = [total_walk.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, ws, dtype=int)[:-1]]
        
        lin_approx = []
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            
            # Get the linear regression parameters
            slope, intercept, _, _, _ = linregress(np.arange(len(ts)), ts)
            
            # Linear approximation
            line = [intercept + slope * dom for dom in np.arange(len(ts))]
            
            # Necessary
            assert len(ts) == len(line)
            
            # Create the linear approximation one
            lin_approx.append(line)
        
        # So we can flatten it
        lin_approx = np.array(lin_approx).flatten()
        partial_time_series = np.array(partial_time_series).flatten()
        
        # Necessary
        assert len(partial_time_series) == len(lin_approx)
        
        plt.plot(partial_time_series)
        plt.plot(lin_approx)
        plt.title(ws)
        plt.show()
        plt.close()
        
        # Mean Squared Deviation of the linear approximation
        f = np.sqrt(((partial_time_series - lin_approx) ** 2).sum() / len(lin_approx))
        
        # Add it to the Mean Squared Error of the window size ws
        mse_array.append(f)
    
    # So that we can regress successfully
    assert len(window_sizes) == len(mse_array)
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(window_sizes), np.log(mse_array))
    
    # See what it looks like
    plt.plot(window_sizes, mse_array)
    plt.plot(window_sizes, [np.exp(intercept) * (l ** slope) for l in window_sizes])
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    plt.close()
    exit()
    
    return [window_sizes, mse_array, slope, intercept, std_err]


def hurst_short(lineage, min_ws, max_ws=None, window_size_steps=5, steps_between_windows=2):
    """ This is following the notation and equations on Wikipedia """
    
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Maximum window size possible
    if not max_ws:
        max_ws = total_length / 2
    
    # What are the window sizes? In order for it not to take a long time we skip two
    window_sizes = np.arange(min_ws, max_ws, window_size_steps, dtype=int)
    
    # The regression of the log of this and the log of the window sizes is the Hurst exponent
    y_to_regress = []
    
    for ws in window_sizes:
        
        # The windowed/partial time-series
        partial_time_series = [lineage.iloc[starting_point:starting_point + ws].values for starting_point in np.arange(0, total_length, steps_between_windows, dtype=int)[:-1]]
        
        # This is where we keep the rescaled ranges for all partial time-series
        rescaled_ranges_array = []
        
        # Go through all the partial time-series
        for ts in partial_time_series:
            # centered time-series
            y = ts - (ts.sum() / len(ts))
            
            # Their walk around the mean
            walk = y.cumsum()
            
            # The absolute range
            ts_range = walk.max() - walk.min()
            
            # For scaling
            ts_std = walk.std()
            
            # Sometimes this happens because inter-division time is discrete and if we have the same value
            # for a certain amount of generations in a row then of course the std is zero.
            if ts_std == 0:
                print('standard div of this partial time-series was zero, so we exclude it from the average')
                continue
            
            # Scaled absolute range
            ts_rescaled_range = ts_range / ts_std
            
            # add it to the array
            rescaled_ranges_array.append(ts_rescaled_range)
            
            # print(y, walk, ts_range, ts_std, ts_rescaled_range, sep='\n')
            # exit()
        
        # Get the average rescaled range
        average_rescaled_range = np.array(rescaled_ranges_array).mean()
        
        # Append it to the array
        y_to_regress.append(average_rescaled_range)
    
    # Convert it to a numpy array
    y_to_regress = np.array(y_to_regress)
    
    # print(window_sizes)
    # print(y_to_regress)
    # print(len(window_sizes))
    # print(len(y_to_regress))
    # exit()
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(window_sizes), np.log(y_to_regress))
    
    return [window_sizes, y_to_regress, slope, intercept, std_err]


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
    seaborn_preamble()
    
    # Histogram of H of every kind
    for kind in histogram_of_regression.kind.unique():
        sns.boxplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=symbols['physical_units'].values(), showfliers=False)
        plt.legend(title='')
        plt.title(kind)
        plt.ylabel('exponent')
        plt.xlabel('')
        plt.savefig('{}/{} per variable.png'.format(figure_folder, kind), dpi=300)
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
        create_folder('{}/Together Scalings'.format(figure_folder))
        for kind in loglog_plot_recreation.kind.unique():
            create_folder('{}/Together Scalings/{}'.format(figure_folder, kind))
    
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
                    (loglog_plot_recreation['kind'] == kind) & (loglog_plot_recreation['variable'] == variable) & (loglog_plot_recreation['dataset'].isin(['Trace', 'Artificial']))].copy()
                
                seaborn_preamble()
                
                set_trace_legend, set_artificial_legend = False, False
                
                wss_array = {'Trace': np.array([]), 'Artificial': np.array([])}
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
                        # The constant artificial color
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
                    # If it is a artificial lineage plot it and include it in the legend
                    elif relevant.loc[ind]['dataset'] == 'Artificial':
                        # The constant artificial color
                        color = cmap[1]
                        
                        # Add the analysis curve to the array of its dataset for the population regression
                        wss_array['Artificial'] = np.append(wss_array['Artificial'], wss)
                        ytr_array['Artificial'] = np.append(ytr_array['Artificial'], ytr)
                        
                        # Include it in the legend but not more than once
                        if set_artificial_legend:
                            plt.plot(wss, ytr, color=color, marker='x')
                            pass
                        else:
                            plt.plot(wss, ytr, color=color, marker='x', label='Artificial')
                            set_artificial_legend = True
                    # We do not want to see the white noise
                    else:
                        continue
                
                # Get the linear regression of all the trajectories pooled for each dataset
                slope_trace, intercept_trace, _, _, std_err_trace = linregress(np.log(np.array(wss_array['Trace']).flatten()), np.log(np.array(ytr_array['Trace']).flatten()))
                slope_art, intercept_art, _, _, std_err_art = linregress(np.log(np.array(wss_array['Artificial']).flatten()), np.log(np.array(ytr_array['Artificial']).flatten()))
                
                # Plot the best fit line and its parameters
                plt.plot(np.unique(wss_array['Trace']), [np.exp(intercept_trace) * (l ** slope_trace) for l in np.unique(wss_array['Trace'])], ls='--', color='blue', linewidth=3,
                         label=r'$' + str(np.round(intercept_trace, 2)) + r'n^{' + str(np.round(slope_trace, 2)) + r'\pm' + str(np.round(std_err_trace, 3)) + r'}$')
                plt.plot(np.unique(wss_array['Artificial']), [np.exp(intercept_art) * (l ** slope_art) for l in np.unique(wss_array['Artificial'])], ls='--', color='red', linewidth=3,
                         label=r'$' + str(np.round(intercept_art, 2)) + r'n^{' + str(np.round(slope_art, 2)) + r'\pm' + str(np.round(std_err_art, 3)) + r'}$')
                
                plt.title(variable)
                plt.ylabel(r'$F(n)$')
                plt.xlabel('n')
                plt.yscale('log')
                plt.xscale('log')
                plt.legend(title='')
                plt.tight_layout()
                # plt.show()
                plt.savefig('{}/Together Scalings/{}/{}.png'.format(figure_folder, kind, variable), dpi=300)
                plt.close()
    # plt.scatter(window_sizes, y_to_regress, color='blue')


def main(args):
    # # the parameters of the scaling analysis
    # if args.data_origin == 'SM':
    #     min_ws = 5
    #     max_ws = None
    #     window_size_steps = 1
    #     steps_between_windows = 1
    # else:
    #     min_ws = 5
    #     max_ws = None
    #     window_size_steps = 2
    #     steps_between_windows = 3

    min_ws = 5
    max_ws = None
    window_size_steps = 2
    steps_between_windows = 3
    
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

            # If we have enough points to get a slope using this lineage
            if len(np.arange(min_ws, len(trace) / 2, window_size_steps, dtype=int)) < 3:
                print('not enough lineages')
                continue
            
            # What type of lineages do we want to look at?
            types_of_lineages, names_types_of_lineages = [trace, artificial, white_noise], ["Trace", "Artificial", "White Noise"]
            
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
        # ax.annotate('Trace: {:.2}\nArtificial: {:.2}'.format(trace_avg_h, artificial_avg_h), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
        #             bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
        
        print('*' * 200)
    
    print('scaling exponents', histogram_of_regression, sep='\n')
    print('recreate df', loglog_plot_recreation, sep='\n')
    
    # Save them
    histogram_of_regression.to_csv('{}/scaling_exponents.csv'.format(args.save_folder), index=False)
    loglog_plot_recreation.to_csv('{}/loglog_scaling_recreation.csv'.format(args.save_folder), index=False)
    
    print('finished saving them')
    
    # # Plot the figures
    # plot_histograms_of_scaling_exponents(args.figs_location, histogram_of_regression)
    # recreate_loglog_plots(loglog_plot_recreation, args.figs_location, individuals=True)


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # # Do all the Mother Machine data
    # for data_origin in mm_data_names:
    #     print(data_origin)
    #
    #     data_origin == 'lambda_LB'
    #
    #     # if data_origin == 'lambda_LB':
    #     #     continue
    #
    #     parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
    #     parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    #     parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
    #                         required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    #     parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
    #                         required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    #     args = parser.parse_args()
    #
    #     main(args)
    #
    #     df = pd.read_csv('{}/scaling_exponents.csv'.format(args.save_folder))
    #
    #     plot_histograms_of_scaling_exponents(args.figs_location, df)
    #
    #     df = pd.read_csv('{}/loglog_scaling_recreation.csv'.format(args.save_folder))
    #
    #     recreate_loglog_plots(df, args.figs_location, individuals=False)
    #
    #     print('*' * 200)
    
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
    
    df = pd.read_csv('{}/scaling_exponents.csv'.format(args.save_folder))
    
    plot_histograms_of_scaling_exponents(args.figs_location, df)
    
    df = pd.read_csv('{}/loglog_scaling_recreation.csv'.format(args.save_folder))
    
    recreate_loglog_plots(df, args.figs_location, individuals=False)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
