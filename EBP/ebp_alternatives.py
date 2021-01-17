#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations

""" vd conditioning on trap as well as lineage """


def vd_with_trap(args):

    pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    tc = pd.read_csv(args['tc']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    
    global_pu = pu[phenotypic_variables].mean()
    
    var_delta = pd.DataFrame(columns=phenotypic_variables)
    var_diff = pd.DataFrame(columns=phenotypic_variables)
    var_tmean = pd.DataFrame(columns=phenotypic_variables)

    array1_delta = pd.DataFrame(columns=phenotypic_variables)
    array1_diff = pd.DataFrame(columns=phenotypic_variables)
    array1_trap = pd.DataFrame(columns=phenotypic_variables)
    
    for trap_id in pu.trap_ID.unique():
        t_cond = (pu['trap_ID'] == trap_id)
        trap = pu[t_cond].copy()
        trap_means = trap[phenotypic_variables].mean()
        
        for trace in ['A', 'B']:
            lin_cond = (trap['trace'] == trace)
            lin = trap[lin_cond].copy()
            
            array1_trap = array1_trap.append(lin.count() * ((trap_means - global_pu) ** 2), ignore_index=True)
            array1_diff = array1_diff.append(lin.count() * ((trap_means - lin[phenotypic_variables].mean()) ** 2), ignore_index=True)
            array1_delta = array1_delta.append(((lin[phenotypic_variables] - lin[phenotypic_variables].mean()) ** 2).sum(), ignore_index=True)
    
    delta_var = array1_delta.sum() / (pu[phenotypic_variables].count() - 1)
    diff_var = array1_diff.sum() / (pu[phenotypic_variables].count() - 1)
    tmean_var = array1_trap.sum() / (pu[phenotypic_variables].count() - 1)
    
    print(delta_var)
    print(tc[phenotypic_variables].var())
    

""" Plot the variance/covariance decompositions and correlations between phenotypic variables as pyramid heatmaps """


def pyramid_heatmaps(args, annot, input_args):
    # pyramid_of_pairwise_covariances()
    def normalize_correlation(df, variables, kind):
        cov = df.cov()
        corr_df = pd.DataFrame(columns=variables, index=variables, dtype=float)
        for param_1 in variables:
            for param_2 in variables:
                if kind == 'decomposition':
                    # It is NOT a pearson correlation, but is instead normalized by the standard deviations of the pooled ensemble
                    normalization = (physical_units[param_1].std() * physical_units[param_2].std())
                elif kind == 'pearson':  # PEARSON CORRELATION
                    normalization = (df.drop_duplicates()[param_1].std() * df.drop_duplicates()[param_2].std())
                else:
                    raise IOError('wrong kind of normalization')
                
                # # If the value lies in the noise range then make it zero
                # if -.1 < cov.loc[param_1, param_2] / normalization < .1:
                #     corr_df.loc[param_1, param_2] = float('nan')  # 0
                # else:
                #     corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
                corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
        
        return corr_df
    
    # read the csv file where we keep the data
    # import/create the trace lineages
    physical_units = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    population_sampled = shuffle_info(physical_units, mm=args['MM'])
    
    for type_of_var, variables in input_args.variable_mapping.items():
        
        mask = np.ones_like(normalize_correlation(physical_units, variables, kind='decomposition'))
        mask[np.tril_indices_from(mask)] = False
        vmax, vmin = 1, -1
        
        level1 = args['ebp_folder'] + '/' + type_of_var
        create_folder(level1)
        
        print('type of var:', type_of_var)
        
        for kind in input_args.kinds_of_correlations:
            
            level2 = level1 + '/' + kind
            create_folder(level2)
            
            print('kind:', kind)
            for label in ['Trace', 'Artificial']:
                if label == 'Trace':
                    pu = physical_units
                    ta = get_time_averages_df(physical_units, phenotypic_variables)
                    tc = trace_center_a_dataframe(physical_units, args['MM'])
                else:
                    pu = population_sampled
                    ta = get_time_averages_df(population_sampled, phenotypic_variables)
                    tc = trace_center_a_dataframe(population_sampled, args['MM'])
                
                print('label:', label)
                
                seaborn_preamble()
                
                if type_of_var == 'phenotypic_variables':
                    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[7 * 2, 2.5 * 2])
                else:
                    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[7 * 1.5, 2.5 * 1.5])
                
                pu = normalize_correlation(pu, variables, kind=kind).rename(columns=symbols['physical_units'], index=symbols['physical_units'])
                ta = normalize_correlation(ta, variables, kind=kind).rename(columns=symbols['time_averages'], index=symbols['time_averages'])
                tc = normalize_correlation(tc, variables, kind=kind).rename(columns=symbols['trace_centered'], index=symbols['trace_centered'])
                
                sns.heatmap(pu, annot=annot, center=0, vmax=vmax,
                            vmin=vmin, cbar=False, ax=axes[0], mask=mask, square=True, fmt='.2f')
                axes[0].set_title('A', x=-.2, fontsize='xx-large')
                
                sns.heatmap(ta, annot=annot, center=0, vmax=vmax, vmin=vmin,
                            cbar=False, ax=axes[1], mask=mask, square=True, fmt='.2f')
                axes[1].set_title('B', x=-.2, fontsize='xx-large')
                
                cbar_ax = fig.add_axes([.91, .1, .03, .8])
                
                sns.heatmap(tc, annot=annot, center=0, vmax=vmax,
                            vmin=vmin, cbar=True, ax=axes[2], mask=mask, cbar_kws={"orientation": "vertical"}, square=True, cbar_ax=cbar_ax, fmt='.2f')
                axes[2].set_title('C', x=-.2, fontsize='xx-large')
                
                plt.suptitle('{} of {} lineage: '.format(kind, label) + args['data_origin'])
                
                fig.tight_layout(rect=[0, 0, .9, 1])
                
                plt.savefig('{}/{}_{}.png'.format(level2, label, args['data_origin']), dpi=300)
                # plt.show()
                plt.close()
                
                print('done!')


""" Shows how the variance decomposition of the lineage decreases per generation """


def ebp_per_gen(df_dict, folder_name):
    seaborn_preamble()
    
    mark = {key: val for key, val in zip(df_dict.keys(), ['o', '^', 'x', '*'])}
    
    repeats = []
    for param1 in phenotypic_variables:
        for param2 in phenotypic_variables:
            if param2 not in repeats:
                for key, df in df_dict.items():
                    print(key)
                    relevant = df[(df['param1'] == param1) & (df['param2'] == param2)]
                    relevant.loc[:, 'generation'] = relevant['generation'] + 1
                    relevant = relevant[relevant['generation'] < 16]
                    
                    # sns.boxplot(data=relevant, x='generation', y='gamma_ta', hue='label', showfliers=False)
                    color_ind = 0
                    slope_array = []
                    for lll in relevant.label.unique():
                        for_reg = relevant[relevant['label'] == lll].copy()
                        # print(for_reg.generation.unique())
                        # print(relevant.generation.unique())
                        # exit()
                        # x_axis = np.append(relevant.generation.unique()[0], relevant.generation.unique() + 1)
                        slope, intercept, _, _, std_err = linregress(np.log(for_reg.generation.values), np.log(for_reg.gamma_ta.values))
                        
                        # sns.lineplot(relevant.generation.unique(), [np.exp(intercept) * ((gen) ** slope) for gen in relevant.generation.unique()], ls='--', color=cmap[color_ind], alpha=.5)
                        color_ind += 1
                        
                        slope_array.append(slope)
                        if np.isnan(slope):
                            
                            slope, intercept, _, _, std_err = linregress(np.log(for_reg.generation.values), np.log(-for_reg.gamma_ta.values))
                            if np.isnan(slope) or np.exp(intercept) < .2:  # It is more or less zero so we don't regress   or np.array([np.isnan(s) for s in slope_array]).any()
                                pass
                            else:  # It was just minus
                                relevant.loc[relevant['label'] == lll, 'label'] = [l + r'$: -{:.2} \, n^'.format(np.exp(intercept)) + '{' + r'{:.2}\pm{:.1}'.format(slope, std_err) + '}$' for l in
                                                                                   for_reg['label'].values]
                                
                                plt.plot(relevant.generation.unique() - 1, [-np.exp(intercept) * ((gen) ** slope) for gen in relevant.generation.unique()], ls='--', alpha=.8)
                        else:
                            relevant.loc[relevant['label'] == lll, 'label'] = [l + r'$: {:.2} \, n^'.format(np.exp(intercept)) + '{' + r'{:.2}\pm{:.1}'.format(slope, std_err) + '}$' for l in
                                                                               for_reg['label'].values]
                            plt.plot(relevant.generation.unique() - 1, [np.exp(intercept) * ((gen) ** slope) for gen in relevant.generation.unique()], ls='--', alpha=.8)
                    
                    sns.pointplot(data=relevant, x='generation', y='gamma_ta', hue='label', marker=mark[key], ci=95, cmap=cmap, label=key, dodge=True,
                                  join=False)
                
                plt.legend(title='')
                # plt.loglog()
                if param1 == param2:
                    plt.ylabel(r'$\Gamma_n(${}$)$'.format(symbols['new_pu'][param1]))
                else:
                    plt.ylabel(r'$\Gamma_n(${}$, ${}$)$'.format(symbols['new_pu'][param1], symbols['new_pu'][param2]))
                plt.axhline(0, color='black')
                plt.yticks(ticks=np.arange(-1.1, 1, .1), labels=['{}'.format(np.round(num, 1)) for num in np.arange(-1.1, 1, .1)])
                plt.ylim([-1, 1])
                plt.xlabel('n')
                # plt.xlim(right=9)
                plt.savefig(folder_name + '/' + param1 + ', ' + param2 + '.png', dpi=300)
                # plt.show()
                plt.close()
        
        repeats.append(param1)


""" Shows the bootstrapped coefficient of variation of all Trace and Population lineages per phenotypic variable. """


def ergodicity_per_variable(eb_df, variable_pairs):
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    _, ax = plt.subplots(tight_layout=True, figsize=[7 * 1.5, 2.5 * 1.5])
    
    necessary_variables = list(np.unique(variable_pairs))
    
    return_proper_latex = lambda param1, param2: symbols['time_averages'][param1] if param1 == param2 else r'({}, {})'.format(symbols['time_averages'][param1], symbols['time_averages'][param2])
    
    the_order = [return_proper_latex(param1, param2) for count, param1 in
                 enumerate(necessary_variables) for param2 in necessary_variables if param2 not in necessary_variables[:count]]
    
    # The latex labels instead of the variable names
    to_plot = eb_df[(eb_df['param1'].isin(necessary_variables)) | (eb_df['param2'].isin(necessary_variables))].copy()
    
    sns.barplot(x='variable', y='gamma_ta', data=to_plot, hue='label', order=the_order, ax=ax, edgecolor='black')  # , label=r'$\overline{cov}$ \ $\sigma^2$'
    ax.yaxis.grid(True)
    ax.set_xlabel('')
    ax.set_ylabel(r'$\Gamma$')
    ax.get_legend().remove()
    plt.show()
    plt.close()


""" Plot the ratio of variance-decompositions/correlations between trace and artificial lineages for each dataset represented as a point """


def ratio_of_all_datasets(datasets, variable_pairs):
    # where we will put all the ratios
    ratio_df = pd.DataFrame(columns=['dataset', 'param1', 'param2', 'variable', 'gamma_ta_ratio', 'pearson_ratio'])
    
    necessary_variables = list(np.unique(variable_pairs))
    
    return_proper_latex = lambda param1, param2: symbols['time_averages'][param1] if param1 == param2 else r'({}, {})'.format(symbols['time_averages'][param1], symbols['time_averages'][param2])
    
    for ds in datasets:
        print(ds)
        df = pd.read_csv(os.path.dirname(os.path.abspath(__file__)) + '/{}/total_length_vd.csv'.format(ds)).drop_duplicates()
        
        # df = df.replace(dict(zip([r'$({}, {})$'.format(symbols['time_averages'][param1], symbols['time_averages'][param2]) for count, param1 in
        #                           enumerate(necessary_variables) for param2 in necessary_variables if param2 not in necessary_variables[:count]],
        #                          [r'('+symbols['time_averages'][param1]+r', ' + symbols['time_averages'][param2] + r')' for count, param1 in
        #                           enumerate(necessary_variables) for param2 in necessary_variables if param2 not in necessary_variables[:count]])))
        
        # print(df)
        # print(df.columns)
        
        for param1, param2, variable in df[['param1', 'param2', 'variable']].drop_duplicates().values:
            rel = df[(df['variable'] == variable)]
            g_ratio = rel[rel['label'] == 'Trace'].gamma_ta.values / rel[rel['label'] == 'Artificial'].gamma_ta.values
            p_ratio = 0  # rel[rel['label'] == 'Trace']['pearson'].values / rel[rel['label'] == 'Artificial']['pearson'].values
            ratio_df = ratio_df.append({
                'dataset': ds, 'param1': param1, 'param2': param2, 'variable': variable, 'gamma_ta_ratio': g_ratio[0], 'pearson_ratio': p_ratio  # p_ratio[0]
            }, ignore_index=True)
    
    # the_order = [return_proper_latex(param1, param2) for count, param1 in
    #              enumerate(necessary_variables) for param2 in necessary_variables if param2 not in necessary_variables[:count]]
    
    # the_order = [r'('+symbols['time_averages'][param1]+r', ' + symbols['time_averages'][param2] + r')' for param1, param2 in variable_pairs]
    the_order = [symbols['time_averages'][param1] for param1, param2 in variable_pairs]
    
    print('the order', the_order)
    print(ratio_df.variable.unique())
    
    # the_order = [return_proper_latex(param1, param1) for count, param1 in enumerate(phenotypic_variables)]
    
    # The latex labels instead of the variable names
    # to_plot = ratio_df[(ratio_df['param1'].isin(necessary_variables)) | (ratio_df['param2'].isin(necessary_variables))].copy()
    to_plot = ratio_df[(ratio_df['variable'].isin(the_order))].copy()
    
    seaborn_preamble()
    
    # sns.barplot(x='variable', y='gamma_ta_ratio', data=to_plot, hue='dataset', order=the_order, edgecolor='black')  # , label=r'$\overline{cov}$ \ $\sigma^2$'
    sns.stripplot(x='variable', y='gamma_ta_ratio', data=to_plot, hue='dataset', order=the_order, linewidth=1, jitter=True)  # , label=r'$\overline{cov}$ \ $\sigma^2$', order=the_order,
    # for ds in datasets:
    #     sns.lineplot(x='variable', y='gamma_ta_ratio', data=to_plot[to_plot['dataset'] == ds], hue='dataset', dodge=True)
    # sns.lineplot(x='variable', y='gamma_ta_ratio', data=to_plot, hue='dataset')
    # sns.pointplot(x='variable', y='gamma_ta_ratio', data=to_plot[['variable', 'gamma_ta_ratio']], order=the_order, color='black', ci='sd', markers='s', capsize=.2)
    # plt.grid(True)
    # plt.yscale('log')
    plt.xlabel('')
    plt.ylabel(r'$\frac{\Gamma_{Trace}}{\Gamma_{Artificial}}$')
    plt.legend('')
    # plt.get_legend().remove()
    # plt.savefig('ratios_of_gamma_ta.png', dpi=300)
    plt.show()
    plt.close()


""" Creates a dataframe with the variance/covariance decomposition of each phenotypic variable """


def create_vd_dataframe(df, phenotypic_variables, kind, n_boots=0):  # MM
    
    # Initialize where we will put the bootstrapped ergodicity breaking variable
    eb_df = pd.DataFrame(columns=['variable', 'kind', 'value'])
    
    # get the dataframes where in each entry, instead of the generation specific value, there is the time-average
    time_averages = get_time_averages_df(df, phenotypic_variables)  # MM
    
    # to normalize the different phenotypic_variables
    pop_var = df.var()
    
    # bootstrap this ergodicity breaking parameter
    if n_boots != 0:
        for _ in np.arange(n_boots):
            # Bootstrapping the indices has the lineage length inequality taken into account
            indices = time_averages.sample(frac=1, replace=True).index
            
            # Get the variance of the time-averages for both kinds of lineages
            variance = time_averages[phenotypic_variables].loc[indices].var() / pop_var
            
            # add them both to the dataframe where we save it all
            for param in phenotypic_variables:
                eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    else:
        # Get the variance of the time-averages for both kinds of lineages
        variance = time_averages[phenotypic_variables].var() / pop_var
        
        # add them both to the dataframe where we save it all
        for param in phenotypic_variables:
            eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    
    return [time_averages, eb_df]


""" adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """


def expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, label, out_df, bs, total_length_df, total_length_bs_df, first_generation, final_generation, bootstraps=False):
    # For the pooled ensemble
    pooled_ensemble_total = pd.DataFrame(columns=df.columns)
    
    # create the info dataframe with values being the cumulative mean of the lineage
    for lin_id in df.lineage_ID.unique():
        # print(lin_id)
        
        # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
        lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
        
        # print(lineage)
        # exit()
        
        # print(len(lineage.generation))
        # print(len(lineage.generation.unique()))
        
        # add its time-average up until and including the total number of generations without NaNs
        to_add = pd.DataFrame.from_dict({
            'label': [label for _ in range(len(lineage))],
            'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
            'generation': [gen for gen in np.arange(len(lineage.generation.unique()))]
        }).reset_index(drop=True)
        
        # Get the expanding mean and the cumulative sum, make sure they are robust to the first NaN in division ratio
        to_add_expanding_mean = pd.concat([lineage[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
        expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
        
        if to_add_expanding_mean.isnull().values.any():
            raise IOError('Got NaNs in the expanding mean, meaning something is wrong')
        
        lineage['generation'] = np.arange(len(lineage), dtype=int)
        
        pooled_ensemble_total = pooled_ensemble_total.append(lineage, ignore_index=True)
        
        if len(to_add_expanding_mean) != len(lineage):
            
            print(lineage)
            print(to_add_expanding_mean)
            raise IOError('Wronfg lenghts')
        
        if not np.array_equal(to_add_expanding_mean['generation'].values, lineage['generation'].values):
            print(lin_id)
            
            print(type(lineage['generation'].iloc[0]))
            print(type(to_add_expanding_mean['generation'].iloc[0]))
            print(lineage['generation'].iloc[0])
            print(to_add_expanding_mean['generation'].iloc[0])
            raise IOError('Wronfg values of generation')
    
    print('out of lineage ID loop')
    
    # This is important
    assert not expanding_mean.isnull().values.any()
    
    # This is important
    assert not pooled_ensemble_total.isnull().values.any()
    
    print('starting time averages')
    
    # Get the EBP for the total length of the lineage
    time_averages = get_time_averages_df(df, phenotypic_variables)
    
    print('time_averages is done')
    
    # for each parameter combination decompose the pooled variance
    tl_temp = {}
    tl_bs_temp = {}
    ind_tl_bs, ind_tl = 0, 0
    repeat = []
    for param1 in phenotypic_variables:
        for param2 in phenotypic_variables:
            if param2 not in repeat:
                ta_cov = time_averages[[param1, param2]].dropna().cov().values[0, 1]
                
                if np.isnan(ta_cov):
                    print(time_averages[time_averages.isnull()])
                    print(')' * 300)
                    print(time_averages[param1].isnull().values.any())
                    print(time_averages[param2].isnull().values.any())
                    print(ta_cov)
                    print(time_averages[param1].values, time_averages[param2].values)
                    print('error with nan')
                    exit()
                
                pool_var = df[param1].std() * df[param2].std()
                
                # Calculate the variance of the time-averages
                gamma_ta_cov = ta_cov / pool_var
                
                # Add it to the dataframe to output
                tl_temp[ind_tl] = {
                    'param1': param1, 'param2': param2,
                    'variable': symbols['time_averages'][param1] if param1 == param2 else r'({}, {})'.format(symbols['time_averages'][param1], symbols['time_averages'][param2]),
                    'n': len(time_averages[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label,
                    'pearson': time_averages[[param1, param2]].dropna().drop_duplicates().corr().values[0, 1]
                }
                
                # print('tl temp dict')
                
                ind_tl += 1
                
                if bootstraps != False:
                    for _ in np.arange(bootstraps):
                        # sample randomly from the distribution with replacement
                        time_averages_replacement = time_averages[[param1, param2]].sample(replace=True, frac=1).dropna().reset_index(drop=True).copy()
                        
                        # Get the covariance of each type of series
                        ta_cov = time_averages_replacement.cov().values[0, 1]  # np.cov(time_averages[param1].values, time_averages[param2].values)[0, 1]
                        
                        pool_var = df[param1].std() * df[param2].std()
                        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                        
                        # Add it
                        tl_bs_temp[ind_tl_bs] = {
                            'param1': param1, 'param2': param2,
                            'variable': symbols['time_averages'][param1] if param1 == param2 else r'({}, {})'.format(symbols['time_averages'][param1], symbols['time_averages'][param2]),
                            'n': len(time_averages_replacement[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label,
                            'pearson': time_averages_replacement.drop_duplicates().corr().values[0, 1]
                        }
                        
                        ind_tl_bs += 1
        
        # To not repeat the calculation twice
        repeat.append(param1)
        
        # print('append param1')
        
        # Append them to the big dataframes
        total_length_df = total_length_df.append(pd.DataFrame.from_dict(tl_temp, "index"), ignore_index=True)
        total_length_bs_df = total_length_bs_df.append(pd.DataFrame.from_dict(tl_bs_temp, "index"), ignore_index=True)
        
        # print('added')
    
    # Important
    if total_length_df.isnull().values.any():
        print('ERROR! total_length_df')
        print(total_length_df)
        exit()
    if total_length_bs_df.isnull().values.any():
        print('ERROR! total_length_bs_df')
        print(total_length_bs_df)
        exit()
    
    # # Calculate the ergodicity breaking parameter over all lineages in the dataset per generation
    # for generation in np.arange(first_generation, final_generation + 1):
    #     print(generation)
    #
    #     out_temp = {}
    #     bs_temp = {}
    #     ind_bs, ind_out = 0, 0
    #
    #     # Get the time-averages of the lineages that have the generation we are looking for
    #     time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)].copy()  # Excludes all lineages that do not reach this generation
    #
    #     # Get the array of which lineages are long enough, ie. have the amount of cycles we need
    #     long_enough_lineages = time_averages['lineage_ID'].unique()
    #
    #     # Define the pooled ensemble of all the lineages that have at least this generation
    #     pooled_ensemble = pooled_ensemble_total[(pooled_ensemble_total['generation'] <= generation) & (pooled_ensemble_total['lineage_ID'].isin(long_enough_lineages))].copy()
    #
    #     # That way we have a robust variance across lineages
    #     if len(long_enough_lineages) > 1:
    #
    #         # for each parameter combination decompose the pooled variance
    #         repeat = []
    #         for param1 in phenotypic_variables:
    #             for param2 in phenotypic_variables:
    #                 if param2 not in repeat:
    #
    #                     # The mean is a linear operator
    #                     if np.abs(pooled_ensemble[param2].mean() - time_averages[param2].mean()) > 0.000001:
    #                         print(pooled_ensemble[param2].mean(), time_averages[param2].mean())
    #                         print(pooled_ensemble[param2].mean() - time_averages[param2].mean())
    #                         raise IOError('Means are not the same in a big way')
    #
    #                     # We must have the same amount of
    #                     if not ((generation + 1) * len(time_averages[param1])) == ((generation + 1) * len(time_averages[param2])) == len(pooled_ensemble[param1]) == len(pooled_ensemble[param2]):
    #                         print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
    #                         raise IOError('Sizes are not the same in a big way')
    #
    #                     # Get the variances of each type of series
    #                     ta_cov = ((generation + 1) * (
    #                             (time_averages[param1].copy() - time_averages.mean()[param1].copy()) * (time_averages[param2].copy() - time_averages.mean()[param2].copy()))).sum() / (
    #                                      len(pooled_ensemble) - 1)
    #                     pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
    #
    #                     # Calculate the variance of the time-averages
    #                     gamma_ta_cov = ta_cov / pool_var
    #
    #                     # Add it to the dataframe to output
    #                     out_temp[ind_out] = {
    #                         'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
    #                     }
    #
    #                     ind_out += 1
    #
    #                     if bootstraps != False:
    #                         for _ in np.arange(bootstraps):
    #                             # sample randomly from the distribution with replacement
    #                             time_averages_replacement = time_averages.sample(replace=True, frac=1)
    #
    #                             # Get the variances of each type of series
    #                             ta_cov = ((generation + 1) * (
    #                                     (time_averages_replacement[param1].copy() - time_averages_replacement.mean()[param1].copy()) * (
    #                                     time_averages_replacement[param2].copy() - time_averages_replacement.mean()[param2].copy()))).sum() / (
    #                                              len(pooled_ensemble) - 1)
    #                             pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
    #
    #                             # Calculate the variance of the time-averages
    #                             gamma_ta_cov = ta_cov / pool_var
    #
    #                             # Add it
    #                             bs_temp[ind_bs] = {
    #                                 'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
    #                             }
    #
    #                             ind_bs += 1
    #
    #             # To not repeat the calculation twice
    #             repeat.append(param1)
    #
    #     # Append them to the big dataframes
    #     out_df = out_df.append(pd.DataFrame.from_dict(out_temp, "index"), ignore_index=True)
    #     bs = bs.append(pd.DataFrame.from_dict(bs_temp, "index"), ignore_index=True)
    
    # Important
    if out_df.isnull().values.any():
        print('ERROR!')
        print(out_df)
        exit()
    if bs.isnull().values.any():
        print('ERROR! BS')
        print(bs)
        exit()
    
    return [out_df, bs, total_length_df, total_length_bs_df]


def main(args, first_generation, final_generation):
    # import/create the trace lineages
    physical_units = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    
    population_sampled = shuffle_info(physical_units, mm=args['MM'])
    
    # lineage_shuffled = shuffle_lineage_generations(physical_units, args['MM'])
    
    # We keep the trap means here
    expanding_mean = pd.DataFrame(columns=['label', 'lineage_ID', 'generation'] + phenotypic_variables)
    
    # Where we keep the gammas
    out_df = pd.DataFrame(columns=['param1', 'param2', 'variable', 'n', 'generation', 'gamma_ta', 'label', 'pearson'])
    bs = pd.DataFrame(columns=['param1', 'param2', 'variable', 'n', 'generation', 'gamma_ta', 'label', 'pearson'])
    total_length_vd = pd.DataFrame(columns=['param1', 'param2', 'variable', 'n', 'generation', 'gamma_ta', 'label', 'pearson'])
    total_length_vd_bs = pd.DataFrame(columns=['param1', 'param2', 'variable', 'n', 'generation', 'gamma_ta', 'label', 'pearson'])
    
    # Calculate the cv and TA per lineage length
    for kind, df in zip(['Trace', 'Artificial'], [physical_units, population_sampled]):
        
        print(kind)
        
        out_df, bs, total_length_vd, total_length_vd_bs = expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, kind, out_df, bs, total_length_vd, total_length_vd_bs,
                                                                                              first_generation, final_generation, bootstraps=False)
        
        # print(out_df, bs, total_length_vd, total_length_vd_bs, sep='\n' * 2)
        
        print('&' * 300)
    
    total_length_vd.to_csv(args['total_length_vd'], index=False)
    # total_length_vd_bs.to_csv('gamma_ta_corrs_bs.csv', index=False)
    out_df.to_csv(args['per_gen_vd'], index=False)
    # bs.to_csv('gamma_ta_corrs_per_gen_bs.csv', index=False)


if __name__ == '__main__':
    import argparse
    import os
    
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
    # parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
    #                     default=dict(zip(['main_variables', 'phenotypic_variables'], [main_variables, phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    kind_of_vd = ['total_length', 'per_gen', 'trap_controlled']
    
    # pairs = list(list(combinations(['division_ratio', 'length_birth', 'generationtime', 'growth_rate'], 2)))  # ['length_birth', 'generationtime', 'growth_rate', 'division_ratio', 'fold_growth']
    # pairs = pairs + [(variable, variable) for variable in list(np.unique(pairs))]  # Add the variances to the covariances : OPTIONAL!
    # pairs = [(variable, variable) for variable in list(np.unique(phenotypic_variables))]
    #
    # ratio_of_all_datasets(sm_datasets[:-1], pairs)
    #
    # exit()
    
    # Do all the Mother and Sister Machine data
    for data_origin in input_args.dataset_names:  # wang_datasets:  # input_args.dataset_names:
        print(data_origin)
        
        create_folder(data_origin)
        
        ebp_folder = os.path.dirname(os.path.abspath(__file__))
        
        processed_data = os.path.dirname(ebp_folder) + '/Datasets/' + data_origin + '/ProcessedData/'
        
        """
        data_origin ==> Name of the dataset we are analysing
        raw_data ==> Where the folder containing the raw data for this dataset is
        processed_data ==> The folder we will put the processed data in
        """
        args = {
            'data_origin': data_origin,
            'raw_data': os.path.dirname(ebp_folder) + '/Datasets/' + data_origin + '/RawData/',
            'MM': False if data_origin in sm_datasets else True,
            # Data singularities, long traces with significant filamentation, sudden drop-offs
            'ebp_folder': ebp_folder,
            'pu': processed_data + 'z_score_under_3/physical_units_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'physical_units.csv',
            'tc': processed_data + 'z_score_under_3/trace_centered_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'trace_centered.csv',
            'total_length_vd': ebp_folder + '/{}/total_length_vd.csv'.format(data_origin),
            'per_gen_vd': ebp_folder + '/{}/per_gen_vd.csv'.format(data_origin),
            'trap_controlled_vd': ebp_folder + '/{}/trap_controlled_vd.csv'.format(data_origin),
            'total_length_figs': ebp_folder + '/{}/total_length_figs/'.format(data_origin),
            'per_gen_figs': ebp_folder + '/{}/per_gen_figs/'.format(data_origin),
            'trap_controlled_figs': ebp_folder + '/{}/trap_controlled_figs/'.format(data_origin),
        }

        vd_with_trap(args)
        exit()
        
        # ebp_with_trap(args)
        #
        # # plot them
        # pyramid_heatmaps(args, annot=True, input_args=input_args)
        
        # '{}variance_decompositions.png'.format(args['total_length_figs'])
        # 
        # copyfile('{}variance_decompositions.png'.format(args['total_length_figs']), 'variance_decompositions_{}.png'.format(args['data_origin']))
        # 
        # continue
        
        #
        # folder_name = args.figs_location + '/ebp_per_gen'
        # create_folder(args.figs_location)
        #
        # main(args, first_generation=0, final_generation=24)
        
        # kinds of vds
        total_length_vd = pd.read_csv(args['total_length_vd'])
        
        pyramid_heatmaps(args, annot=True, input_args=input_args)
        
        #
        pairs = list(list(combinations(['length_birth', 'generationtime', 'growth_rate'], 2)))  # 'fold_growth', 'division_ratio',
        pairs = pairs + [(variable, variable) for variable in list(np.unique(pairs))]  # Add the variances to the covariances : OPTIONAL!
        print(pairs)
        ergodicity_per_variable(total_length_vd, pairs)
        
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        _, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
        sns.barplot(data=total_length_vd, x='variable', y='gamma_ta', hue='label', order=list(symbols['time_averages'].values()), ax=ax,
                    edgecolor='black')  # , label=r'$\overline{cov}$ \ $\sigma^2$'
        ax.yaxis.grid(True)
        ax.set_xlabel('')
        ax.set_ylabel(r'$\Gamma$')
        plt.title(data_origin)
        ax.get_legend().remove()
        # save the figure
        # plt.savefig('{}variance_decompositions.png'.format(args['total_length_figs']), dpi=300)
        plt.show()
        plt.close()
        
        print('*' * 200)
