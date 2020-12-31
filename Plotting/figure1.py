#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units, dataset_names, create_folder
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


""" This creates the illustration of two lineage distributions vs a gray background Population.  """


def kldiv_illustration(info, phenotypic_variables_new, symbols, units, kl_df, type_of_lineage, MM, ax):
    # The latex labels instead of the variable names
    kl_df = kl_df.replace(symbols)
    
    # This is a manual thing to make the illustration look better
    bounds = {
        'generationtime': [.22, 1], 'fold_growth': [.27, 1], 'growth_rate': [.8, 1.8], 'length_birth': [.7, 4.5], 'length_final': [2.3, 8.3], 'division_ratio': [.35, .65], 'added_length': [.4, 5]
    }
    
    # So the colors of the illustrations are not mixed with those of the trace and population
    cmap = sns.color_palette("tab10")[5:7]
    
    for param in phenotypic_variables_new:
        
        # Define the population distribution for the grey background
        pop = info[param]
        
        # Find out which traps have the highest kl_divergences
        sorted_empirical = kl_df[(kl_df['variable'] == symbols[param]) & (kl_df['kind'] == type_of_lineage)].sort_values('value', ascending=False)[['value', 'lineage_ID']]
        
        # Because the generationtime gives us weird binning...
        if param == 'generationtime':
            # Plot the Population distribution, we cut it to get appropriate binning... Ask Naama if ok...
            sns.distplot(pop, bins=12, kde=False, color='gray', norm_hist=True,
                         hist_kws={'range': tuple(bounds[param]), 'alpha': .5, "edgecolor": "white"}, ax=ax)
        else:
            # Plot the Population distribution, we cut it to get appropriate binning... Ask Naama if ok...
            sns.distplot(pop, kde=True, color='gray', norm_hist=True,
                         hist_kws={'range': tuple(bounds[param]), 'alpha': .5, "edgecolor": "white"}, ax=ax)

        # Choose automatically opposing distributions with min and max TA to make see the difference, and make sure it is bigger than 15 generations
        optimal = pd.DataFrame({
            'index': np.arange(0, 25),
            'mean': [info[(info['lineage_ID'] == sorted_empirical.iloc[index]['lineage_ID'])][param].mean() for index in np.arange(0, 25)],
            'size': [len(info[(info['lineage_ID'] == sorted_empirical.iloc[index]['lineage_ID'])][param]) for index in np.arange(0, 25)]
        }).sort_values('mean', ascending=False)
        
        # if MM:
        #     # Choose automatically opposing distributions with min and max TA to make see the difference, and make sure it is bigger than 15 generations
        #     optimal = pd.DataFrame({
        #         'index': np.arange(0, 25),
        #         'mean': [info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID'])][param].mean() for index in np.arange(0, 25)],
        #         'size': [len(info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID'])][param]) for index in np.arange(0, 25)]
        #     }).sort_values('mean', ascending=False)
        # else:
        #     # Choose automatically opposing distributions with min and max TA to make see the difference, and make sure it is bigger than 15 generations
        #     optimal = pd.DataFrame({
        #         'index': np.arange(0, 25),
        #         'mean': [info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[index]['trace'])][param].mean() for index in np.arange(0, 25)],
        #         'size': [len(info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[index]['trace'])][param]) for index in np.arange(0, 25)]
        #     }).sort_values('mean', ascending=False)
        
        # Choose the indices in the kl_div df to use for the illustration
        low, high = optimal[optimal['size'] >= 15].iloc[0]['index'], optimal[optimal['size'] >= 15].iloc[-1]['index']
        
        # For the Reference in the Latex
        lineage_ids = []
        lineage_lengths = []
        
        # Plot the seperate two distributions
        for count, index in enumerate([low, high]):
            # if MM:
            #     # Get the empirical data
            #     empirical = info[(info['trap_ID'] == sorted_empirical.iloc[int(index)]['lineage_ID'])][param]
            # else:
            #     # Get the empirical data
            #     empirical = info[(info['trap_ID'] == sorted_empirical.iloc[int(index)]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[int(index)]['trace'])][param]

            # Get the empirical data
            empirical = info[(info['lineage_ID'] == sorted_empirical.iloc[int(index)]['lineage_ID'])][param]
            
            # For the Reference in the Latex
            lineage_ids.append(sorted_empirical.iloc[int(index)]['lineage_ID'])
            lineage_lengths.append(len(empirical))
            
            # Plot the low empirical
            sns.distplot(empirical, kde=False, norm_hist=True, ax=ax, hist_kws={"edgecolor": "white"}, color=cmap[count])

        ax.set_xlim(bounds[param])
        ax.set_ylim([0, 6])
        ax.set_xlabel(r'{}'.format(symbols[param]) + r' {}'.format(units[param]))
        ax.grid(False)
        ax.set_ylabel('')


""" Shows the distirbutions of the Kl-divergences of all Trace and Population lineages per phenotypic variable. """


def kl_div_per_variable(kl_df, ax, symbols):
    # The latex labels instead of the variable names
    kl_df = kl_df.replace(symbols)
    kl_df = kl_df.replace({'Population': 'Artificial'})
    
    # Here we start plotting both KL-Divergences, The good part with this is that it takes into consideration the standard deviation as well
    # sns.barplot(data=kl_df, x='variable', y='value', hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ci=95, ax=ax)
    # sns.pointplot(data=kl_df, x='variable', y='value', hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'],  ax=ax, join=False, ci=95, dodge=True)
    sns.boxplot(data=kl_df, x='variable', y='value', hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], showfliers=False, ax=ax)
    ax.yaxis.grid(True)
    ax.set_xlabel('')
    ax.set_ylabel(r'$D_{KL}$')
    ax.legend(title='')
    # ax.get_legend().remove()
    ax.set_ylim(bottom=-.01)


""" Shows the bootstrapped coefficient of variation of all Trace and Population lineages per phenotypic variable. """


def ergodicity_per_variable(eb_df, ax, symbolss):
    # The latex labels instead of the variable names
    eb_df = eb_df.replace(symbolss)#.replace({'Trace': , 'Population': })
    
    # total = pd.DataFrame(columns=['variable', 'kind', 'value'])
    # for kind in ['Trace', 'Population']:
    #     for param in [r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$']:
    #         total = total.append({'variable': param, 'value': 1, 'kind': kind}, ignore_index=True)
    
    # plot the cv
    # sns.barplot(x='variable', y='value', data=total, hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ax=ax, edgecolor='black', alpha=.7, hatch = '/') # , label=r'$\delta cov$ \ $\sigma^2$'
    sns.barplot(x='variable', y='value', data=eb_df, hue='kind', order=list(symbolss.values()), ax=ax, edgecolor='black') # , label=r'$\overline{cov}$ \ $\sigma^2$'
    # sns.boxplot(x='variable', y='value', data=eb_df, showfliers=False, hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ax=ax)
    ax.yaxis.grid(True)
    ax.set_xlabel('')
    ax.set_ylabel(r'$\Gamma$')
    ax.get_legend().remove()
    
    
""" This plots the main big figure """
    
    
def main(args):
    # import the labeled measured bacteria in physical units
    info = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    # import the data
    kl_df = pd.read_csv('{}/{}'.format(args.save_folder, args.kld))
    eb_df = pd.read_csv('{}/{}'.format(args.save_folder, args.ebp))
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    
    # create the figure and construct the layout of the figure
    fig = plt.figure(constrained_layout=True, frameon=False, figsize=[7, 7])  # , figsize=[7.5*2, 2.44*2], tight_layout=True
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 0:3])
    ax5 = fig.add_subplot(gs[1, 3:6])
    
    # axis 1: function for the first axes is to delete the frame and ticks of the graph for the illustration
    ax1.set_title('A', x=-.23, fontsize='xx-large')
    ax1.set_frame_on(False)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # axis 2
    kldiv_illustration(info, ['growth_rate'], symbols['new_pu'], units, kl_df, 'Trace', True, ax2)
    
    # axis 3
    kldiv_illustration(info, ['fold_growth'], symbols['new_pu'], units, kl_df, 'Trace', True, ax3)
    
    # axis 4
    kl_div_per_variable(kl_df, ax4, symbols['new_pu'])
    
    # axis 5
    ergodicity_per_variable(eb_df, ax5, symbols['time_averages'])
    
    # Let's put the titles
    ax2.set_title('B', x=-.15, fontsize='xx-large')
    ax4.set_title('C', x=-.15, fontsize='xx-large')
    ax5.set_title('D', x=-.15, fontsize='xx-large')
    
    # save the figure
    plt.savefig('{}/Figure 1.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()
    

def individuals(args):
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
    ergodicity_per_variable(pd.read_csv('{}/{}'.format(args.save_folder, args.ebp)), ax, symbols['time_averages'])
    plt.title(args.data_origin)
    plt.ylim([0, 1])
    # save the figure
    plt.savefig('{}/EBP.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
    kl_div_per_variable(pd.read_csv('{}/{}'.format(args.save_folder, args.kld)), ax, symbols['new_pu'])
    plt.title(args.data_origin)
    # save the figure
    plt.savefig('{}/KLD.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()
    
    
if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # # Do all the Mother Machine data
    # for data_origin in dataset_names:
    #
    #     parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
    #     parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    #     parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
    #                         required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    #     parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
    #                         required=False, default='physical_units_with_outliers.csv')
    #     parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
    #                         required=False, default='population_lineages_with_outliers.csv')
    #     parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
    #                         required=False, default='time_averages_with_outliers.csv')
    #     parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
    #                         required=False, default='ergodicity_breaking_parameter_with_outliers.csv')
    #     parser.add_argument('-kld', '--kld', metavar='', type=str,
    #                         help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
    #                         required=False, default='kullback_leibler_divergences_with_outliers.csv')
    #     parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
    #     parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
    #                         required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    #     args = parser.parse_args()
    #
    #     create_folder(args.figs_location)
    #
    #     # main(args)
    #
    #     # print(pd.read_csv('{}/{}'.format(args.save_folder, args.ebp)))
    #     # exit()
    #
    #     # set a style on seaborn module
    #     sns.set_context('paper')
    #     sns.set_style("ticks", {'axes.grid': True})
    #     _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
    #     ergodicity_per_variable(pd.read_csv('{}/{}'.format(args.save_folder, args.ebp)), ax, symbols['new_pu'])
    #     # save the figure
    #     plt.savefig('{}/EBP.png'.format(args.figs_location), dpi=300)
    #     # plt.show()
    #     plt.close()
    #
    #     # set a style on seaborn module
    #     sns.set_context('paper')
    #     sns.set_style("ticks", {'axes.grid': True})
    #     _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
    #     kl_div_per_variable(pd.read_csv('{}/{}'.format(args.save_folder, args.kld)), ax, symbols['new_pu'])
    #     # save the figure
    #     plt.savefig('{}/KLD.png'.format(args.figs_location), dpi=300)
    #     # plt.show()
    #     plt.close()
    #
    #     print('*' * 200)
    
    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time
    
    data_origin = 'SM'
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units_with_outliers.csv')
    parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                        required=False, default='population_lineages_with_outliers.csv')
    parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                        required=False, default='time_averages_with_outliers.csv')
    parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                        required=False, default='ergodicity_breaking_parameter_with_outliers.csv')
    parser.add_argument('-kld', '--kld', metavar='', type=str,
                        help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                        required=False, default='kullback_leibler_divergences_with_outliers.csv')
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    args = parser.parse_args()
    
    create_folder(args.figs_location)
    
    for blah in ['A', 'B']:
        
        ebp = pd.read_csv('{}/{}_{}'.format(args.save_folder, args.ebp, blah))
    
        # set a style on seaborn module
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
        ergodicity_per_variable(ebp, ax, symbols['time_averages'])
        # save the figure
        plt.ylim([0, .45])
        plt.savefig('{}/EBP_{}.png'.format(args.figs_location, blah), dpi=300)
        # plt.show()
        plt.close()
    
    # main(args)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
