#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, symbols
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
import pingouin as pg


# The function that will do the analysis
def get_deviation_dataframes(variables, out_df, regs, dataframe, label, distances=np.arange(1, 51)):
    for count, lin_id in enumerate(dataframe.lineage_ID.unique()):
        print(count)
        
        # Temporary dictionary to make the process of saving data faster
        temp_dict = {}
        ind = 0
        
        # Get the lineage from the dataframe
        lineage = dataframe[dataframe['lineage_ID'] == lin_id].copy().sort_values('generation')
        
        # Get the cumulative sum for the variables in the lineage
        stand_cumsum_lineage = ((lineage[variables] - lineage[variables].mean()) / lineage[variables].std()).cumsum()
        
        for variable in variables:
            
            # the random walk
            lin_walk = stand_cumsum_lineage[variable].values - stand_cumsum_lineage[variable].dropna().iloc[0]
            
            msd_lineage = []
            
            for distance in distances:
                
                msd = np.mean(np.array([
                    # The calculation
                    (first - second) ** 2
                    
                    # the points in the time-series
                    for first, second in zip(lin_walk[:-distance], lin_walk[distance:])
                    
                    # Only include the points that are not NaNs
                    if ~np.any([np.isnan(first), np.isnan(second)])
                ]))
                
                msd_lineage.append(msd)
                
                temp_dict[ind] = {
                    'MSD': msd,
                    'RMSD': np.sqrt(msd),
                    'NRMSD': np.sqrt(msd) / (np.nanmax(lin_walk) - np.nanmin(lin_walk)),
                    'CV_NRMSD': np.sqrt(msd) / np.nanmean(lin_walk),
                    'distance': distance,
                    'dataset': args.data_origin,
                    'n': len(lin_walk) - distance,
                    'variable': variable,
                    'lineage': lin_id,
                    'kind': label
                }
                
                # make space in the dictionary
                ind += 1
            
            slope, intercept, _, _, std_err = linregress(np.log(distances), np.log(msd_lineage))
            
            regs = regs.append({
                'intercept': np.exp(intercept),
                'slope': slope,
                'std_error': std_err,
                'variable': variable,
                'lineage': lin_id,
                'kind': label
            }, ignore_index=True)
            
            # # Plot it
            # plt.plot(lin_walk)
            # plt.axhline(0, c='black')
            # plt.show()
            # plt.close()
        # print(temp_dict)
        
        # Append them to the big dataframes
        out_df = out_df.append(pd.DataFrame.from_dict(temp_dict, "index"), ignore_index=True)
    
    return [out_df, regs]


def main(args, variables, distances=np.arange(1, 51)):
    # The types of lineages
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    
    population_sampled = shuffle_info(physical_units, mm=args.MM)
    
    lineage_shuffled = shuffle_lineage_generations(physical_units, args.MM)
    
    # This is where we will store the results of our analysis
    mean_square_deviations = pd.DataFrame(columns=['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD', 'distance', 'dataset', 'n', 'variable', 'lineage', 'kind'])
    
    regressions = pd.DataFrame(columns=['intercept', 'slope', 'std_error', 'variable', 'lineage', 'kind'])
    
    # For each type of lineage, do the analysis
    for dataframe, label in zip([physical_units, population_sampled, lineage_shuffled], ['Trace', 'Artificial', 'Shuffled']):
        print(label)
        
        mean_square_deviations, regressions = get_deviation_dataframes(variables, mean_square_deviations, regressions, dataframe, label)
    
    # Save them to the folder
    mean_square_deviations.to_csv('{}/MSD_df.csv'.format(args.data_origin), index=False)
    regressions.to_csv('{}/regression_df.csv'.format(args.data_origin), index=False)
    
    # Make sure all the regressions went well
    assert ~regressions.isnull().values.any()
    
    return [mean_square_deviations, regressions]


def plot_individual_msds(args, variables, dataframe, label):
    create_folder(args.data_origin)
    for y in ['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD']:
        create_folder(args.data_origin + '/' + y)
        create_folder(args.data_origin + '/' + y + '/' + label)
    
    print('plotting')
    
    seaborn_preamble()
    
    for variable in variables:
        # if variable != 'fold_growth':
        #     continue
        
        for y in ['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD']:
            # if y != 'MSD':
            #     continue
            
            relevant = dataframe[(dataframe['variable'] == variable)].sort_values(['lineage', 'distance']).copy()
            
            # slope, intercept, _, _, std_err = linregress(np.log(relevant['distance'].values), np.log(relevant[y].values))
            
            for lin_id in relevant.lineage.unique():
                to_plot = relevant[relevant['lineage'] == lin_id].sort_values('distance')
                # print(to_plot)
                # print(to_plot['distance'])
                # exit()
                plt.plot(to_plot['distance'].values, to_plot[y].values)  # , marker='o' to_plot['distance'].values
            
            # sns.pointplot(data=relevant, x='distance', y='MSD', color='black', ci=95)
            #
            # relevant['distance'] = relevant['distance'] - 1
            # sns.lineplot(data=relevant, x='distance', y=y, ci=None, marker=None, hue='lineage')
            
            the_y_mean = np.array([relevant[relevant['distance'] == dist][y].mean() for dist in relevant.distance.unique()])
            yerr = np.array([relevant[relevant['distance'] == dist][y].std() for dist in relevant.distance.unique()])
            
            plt.errorbar(x=relevant.distance.unique(), y=the_y_mean, yerr=yerr, color='black', marker='o', capsize=1, capthick=1)
            
            # plt.plot(relevant.distance.unique(), [np.exp(intercept) * (distance ** slope) for distance in relevant.distance.unique()], ls='--', linewidth=3,
            #          label=str(np.exp(intercept))[:4] + r'$x^{' + str(slope)[:4] + '}$', color='black')
            
            plt.title(args.data_origin + ': ' + label)
            # plt.legend(title='')
            plt.ylabel(r'{}('.format(y) + symbols['physical_units'][variable] + r')')
            plt.xlabel('Distance')
            # plt.yscale('log')
            # plt.xscale('log')
            plt.tight_layout()
            plt.savefig('{}/{}/{}/{}_nologlog'.format(args.data_origin, y, label, variable), dpi=300)
            # plt.show()
            plt.close()
            


def heatmap_analogs(args, df, label, variables=phenotypic_variables[3:], suffix=''):
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    
    latex_symbols = {variable: symbols[label][variable] for variable in variables}
    
    df = df[variables].copy().rename(columns=latex_symbols)
    
    to_plot = pg.pairwise_corr(df, method='pearson')
    
    # reformat the results into a dataframe we can use for the heatmap
    to_plot_df = pd.DataFrame(columns=to_plot.X.unique(), index=to_plot.Y.unique(), dtype=float)
    repeats = []
    for x in to_plot.X.unique():
        for y in to_plot.Y.unique():
            if x != y and y not in repeats:
                # Add the results
                to_plot_df[x].loc[y] = to_plot[(to_plot['X'] == x) & (to_plot['Y'] == y)]['r'].values[0]
        
        repeats.append(x)
    
    # The mask to show only the lower triangle in the heatmap
    mask = np.ones_like(to_plot_df)
    mask[np.tril_indices_from(mask)] = False
    
    fig, ax = plt.subplots(tight_layout=True, figsize=[7 * (2 / 3), 5.5 * (2 / 3)])
    
    sns.heatmap(to_plot_df, annot=True, square=True, vmin=-1, vmax=1, mask=mask, center=0)
    # plt.title(label)
    # plt.tight_layout()
    # plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scch, label, suffix), dpi=300)
    plt.show()
    plt.close()


def hists_of_regression_stats(args, dataframe):
    name = '{}/Histogram of regression parameters'.format(args.data_origin)
    create_folder(name)
    seaborn_preamble()
    df = dataframe.replace(symbols['physical_units']).copy()
    for y in ['slope', 'intercept', 'std_error']:
        sns.boxplot(data=df, x='variable', y=y, showfliers=False, hue='kind')
        plt.title('')
        plt.xlabel('')
        plt.legend(title='')
        plt.tight_layout()
        # plt.show()
        plt.savefig(name + '/' + y + '.png', dpi=300)
        plt.close()


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # The final_generations for the different datasets
    fg = {
        'SM': 45,
        # 'lambda_LB': 60,
        'Maryam_LongTraces': 45,
        'MG1655_inLB_LongTraces': 200,
        # 'LAC_M9': np.nan
        '4-28-2017': 30
    }
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in ['MC4100_25C', 'MC4100_27C', 'MC4100_37C']:  # dataset_names
        # data_origin = 'lambda_LB'
        print(data_origin)
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                            required=False, default='artificial_lineages.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                            required=False, default='gamma_ta_corrs_per_gen.csv')
        parser.add_argument('-kld', '--kld', metavar='', type=str,
                            help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                            required=False, default='kullback_leibler_divergences.csv')
        parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        create_folder(args.figs_location)
        
        # if data_origin == 'SM':
        #     main(args, phenotypic_variables)  # fg[data_origin]
        # else:
        #     main(args, phenotypic_variables)
            
        print('Finished analysis')
        
        out_df = pd.read_csv('{}/MSD_df.csv'.format(args.data_origin))
        reg_df = pd.read_csv('{}/regression_df.csv'.format(args.data_origin))
        
        # # Plot the boxplots of the regression parameters
        # hists_of_regression_stats(args, reg_df)
        
        # Plot the individual msds for different lineage types and variables
        for label in out_df.kind.unique():
            print(label)
            dataframe = out_df[out_df['kind'] == label].copy()
            
            # Plot it
            plot_individual_msds(args, phenotypic_variables, dataframe, label)
        
        exit()
        # Where to save the correlations between regression parameters
        correlations = pd.DataFrame(columns=['var1', 'var2', 'pcorr', 'reg_prop1', 'reg_prop2', 'kind'])
        
        # Get the correlations between different regression parameters for different lineage types and variables
        for label in out_df.kind.unique():
            print(label)
            
            dataframe = reg_df[reg_df['kind'] == label].copy()
            
            print(dataframe)
            exit()
            
            repeat = []
            for var1 in phenotypic_variables:
                for var2 in phenotypic_variables:
                    if var2 not in repeat:
                        relevant_x = dataframe[dataframe['variable'] == var1]
                        relevant_y = dataframe[dataframe['variable'] == var2]
                        for reg_prop1 in ['slope', 'intercept', 'std_error']:
                            for reg_prop2 in ['slope', 'intercept', 'std_error']:
                                # For the correlation
                                x = relevant_x[reg_prop1].sort_values('lineage').values
                                y = relevant_y[reg_prop2].sort_values('lineage').values
                                
                                correlations.append({
                                    'var1': var1,
                                    'var2': var2,
                                    'pcorr': pearsonr(x, y),
                                    'reg_prop1': reg_prop1,
                                    'reg_prop2': reg_prop2,
                                    'kind': label
                                }, ignore_index=True)
                                print(correlations)
                                exit()
                repeat.append(var1)
        
        print('*' * 200)
        exit()
