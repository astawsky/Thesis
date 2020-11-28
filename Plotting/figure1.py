#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


""" This creates the illustration of two lineage distributions vs a gray background Population.  """


def kldiv_illustration(info, phenotypic_variables, symbols, units, kl_df, type_of_lineage, ax):
    # The latex labels instead of the variable names
    kl_df = kl_df.replace(symbols)
    
    # This is a manual thing to make the illustration look better
    bounds = {
        'generationtime': [.22, 1], 'fold_growth': [.27, 1], 'growth_rate': [.8, 1.8], 'length_birth': [.7, 4.5], 'length_final': [2.3, 8.3], 'division_ratio': [.35, .65], 'added_length': [.4, 5]
    }
    
    # So the colors of the illustrations are not mixed with those of the trace and population
    cmap = sns.color_palette("tab10")[5:7]
    
    for param in phenotypic_variables:
        
        # Define the population distribution for the grey background
        pop = info[param]
        
        # Find out which traps have the highest kl_divergences
        sorted_empirical = kl_df[(kl_df['variable'] == symbols[param]) & (kl_df['kind'] == type_of_lineage)].sort_values('value', ascending=False)[['value', 'trap_ID', 'trace']]
        
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
            'mean': [info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[index]['trace'])][param].mean() for index in np.arange(0, 25)],
            'size': [len(info[(info['trap_ID'] == sorted_empirical.iloc[index]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[index]['trace'])][param]) for index in np.arange(0, 25)]
        }).sort_values('mean', ascending=False)
        
        # Choose the indices in the kl_div df to use for the illustration
        low, high = optimal[optimal['size'] >= 15].iloc[0]['index'], optimal[optimal['size'] >= 15].iloc[-1]['index']
        
        # For the Reference in the Latex
        lineage_ids = []
        lineage_lengths = []
        
        # Plot the seperate two distributions
        for count, index in enumerate([low, high]):
            # Get the empirical data
            empirical = info[(info['trap_ID'] == sorted_empirical.iloc[int(index)]['trap_ID']) & (info['trace'] == sorted_empirical.iloc[int(index)]['trace'])][param]
            
            # For the Reference in the Latex
            lineage_ids.append(sorted_empirical.iloc[int(index)]['trap_ID'])
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


def ergodicity_per_variable(eb_df, ax, symbols):
    # The latex labels instead of the variable names
    eb_df = eb_df.replace(symbols)#.replace({'Trace': , 'Population': })
    
    # total = pd.DataFrame(columns=['variable', 'kind', 'value'])
    # for kind in ['Trace', 'Population']:
    #     for param in [r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$']:
    #         total = total.append({'variable': param, 'value': 1, 'kind': kind}, ignore_index=True)
    
    # plot the cv
    # sns.barplot(x='variable', y='value', data=total, hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ax=ax, edgecolor='black', alpha=.7, hatch = '/') # , label=r'$\delta cov$ \ $\sigma^2$'
    sns.barplot(x='variable', y='value', data=eb_df, hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ax=ax, edgecolor='black') # , label=r'$\overline{cov}$ \ $\sigma^2$'
    # sns.boxplot(x='variable', y='value', data=eb_df, showfliers=False, hue='kind', order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ax=ax)
    ax.yaxis.grid(True)
    ax.set_xlabel('')
    ax.set_ylabel('EB')
    ax.get_legend().remove()
    
    
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
    kldiv_illustration(info, ['growth_rate'], symbols['physical_units'], units, kl_df, 'Trace', ax2)
    
    # axis 3
    kldiv_illustration(info, ['fold_growth'], symbols['physical_units'], units, kl_df, 'Trace', ax3)
    
    # axis 4
    kl_div_per_variable(kl_df, ax4, symbols['physical_units'])
    
    # axis 5
    ergodicity_per_variable(eb_df, ax5, symbols['physical_units'])
    
    # Let's put the titles
    ax2.set_title('B', x=-.15, fontsize='xx-large')
    ax4.set_title('C', x=-.15, fontsize='xx-large')
    ax5.set_title('D', x=-.15, fontsize='xx-large')
    
    # save the figure
    plt.savefig('{}/Figure 1.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()
