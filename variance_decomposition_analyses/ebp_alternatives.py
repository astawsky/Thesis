#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress

""" Returns dataframe with same size containing the time-averages of each phenotypic variable instead of the local value """


def get_time_averages_df(info, phenotypic_variables):  # MM
    # We keep the trap means here
    means_df = pd.DataFrame(columns=['lineage_ID', 'max_gen', 'generation'] + phenotypic_variables)
    
    # specify a lineage
    for trap_id in info['lineage_ID'].unique():
        
        # the values of the lineage we get from physical units
        lineage = info[(info['lineage_ID'] == trap_id)].copy()
        
        # add its time-average
        to_add = {
            'lineage_ID': [trap_id for _ in np.arange(len(lineage))],
            'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
        }
        to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables})
        to_add = pd.DataFrame(to_add)
        means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    
    assert len(info) == len(means_df)
    
    return means_df


""" Two ways to calculate the ergodicity breaking parameter, with an expanding window or using the window of the lineages aleady given """


def get_both_per_gen_and_total_ebp(df, phenotypic_variables, expanding_mean, label, per_gen_ebp, per_gen_ebp_bs, total_length_df, total_length_bs_df, first_generation, final_generation,
                                   bootstraps=False):
    # For the pooled ensemble
    pooled_ensemble_total = pd.DataFrame(columns=df.columns)
    
    # create the info dataframe with values being the cumulative mean of the lineage
    for lin_id in df.lineage_ID.unique():
        
        # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
        lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
        
        # add its time-average up until and including the total number of generations without NaNs
        to_add = pd.DataFrame.from_dict({
            'label': [label for _ in range(len(lineage))],
            'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
            'generation': [generation for generation in np.arange(len(lineage), dtype=int)]
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
            raise IOError('Wronfg values of generation')
    
    # This is important
    assert not expanding_mean.isnull().values.any()
    
    # This is important
    assert not pooled_ensemble_total.isnull().values.any()
    
    # Get the EBP for the total length of the lineage
    time_averages = get_time_averages_df(df, phenotypic_variables)
    
    # for each parameter combination decompose the pooled variance
    tl_temp = {}
    tl_bs_temp = {}
    ind_tl_bs, ind_tl = 0, 0
    repeat = []
    for param1 in phenotypic_variables:
        for param2 in phenotypic_variables:
            if param2 not in repeat:
                # Get the variances of each type of series
                ta_cov = np.cov(time_averages[param1], time_averages[param2])[0, 1]
                
                pool_var = df[param1].std() * df[param2].std()
                
                # Calculate the variance of the time-averages
                gamma_ta_cov = ta_cov / pool_var
                
                # Add it to the dataframe to output
                tl_temp[ind_tl] = {
                    'param1': param1, 'param2': param2, 'n': len(time_averages[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label
                }
                
                ind_tl += 1
                
                if bootstraps != False:
                    for _ in np.arange(bootstraps):
                        # sample randomly from the distribution with replacement
                        time_averages_replacement = time_averages.sample(replace=True, frac=1).copy()
                        
                        # Get the covariance of each type of series
                        ta_cov = np.cov(time_averages[param1], time_averages[param2])[0, 1]
                        
                        pool_var = df[param1].std() * df[param2].std()
                        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                        
                        # Add it
                        tl_bs_temp[ind_tl_bs] = {
                            'param1': param1, 'param2': param2, 'n': len(time_averages[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label
                        }
                        
                        ind_tl_bs += 1
        
        # To not repeat the calculation twice
        repeat.append(param1)
        
        # Append them to the big dataframes
        total_length_df = total_length_df.append(pd.DataFrame.from_dict(tl_temp, "index"), ignore_index=True)
        total_length_bs_df = total_length_bs_df.append(pd.DataFrame.from_dict(tl_bs_temp, "index"), ignore_index=True)
    
    # Important
    if total_length_df.isnull().values.any():
        print('ERROR! total_length_df')
        print(total_length_df)
        exit()
    if total_length_bs_df.isnull().values.any():
        print('ERROR! total_length_bs_df')
        print(total_length_bs_df)
        exit()
    
    # Calculate the ergodicity breaking parameter over all lineages in the dataset per generation
    for generation in np.arange(first_generation, final_generation + 1):
        print(generation)
        
        out_temp = {}
        bs_temp = {}
        ind_bs, ind_out = 0, 0
        
        # Get the time-averages of the lineages that have the generation we are looking for
        time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)].copy()  # Excludes all lineages that do not reach this generation
        
        # Get the array of which lineages are long enough, ie. have the amount of cycles we need
        long_enough_lineages = time_averages['lineage_ID'].unique()
        
        # Define the pooled ensemble of all the lineages that have at least this generation
        pooled_ensemble = pooled_ensemble_total[(pooled_ensemble_total['generation'] <= generation) & (pooled_ensemble_total['lineage_ID'].isin(long_enough_lineages))].copy()
        
        # That way we have a robust variance across lineages
        if len(long_enough_lineages) > 1:
            
            # for each parameter combination decompose the pooled variance
            repeat = []
            for param1 in phenotypic_variables:
                for param2 in phenotypic_variables:
                    if param2 not in repeat:
                        
                        # The mean is a linear operator
                        if np.abs(pooled_ensemble[param2].mean() - time_averages[param2].mean()) > 0.000001:
                            print(pooled_ensemble[param2].mean(), time_averages[param2].mean())
                            print(pooled_ensemble[param2].mean() - time_averages[param2].mean())
                            raise IOError('Means are not the same in a big way')
                        
                        # We must have the same amount of
                        if not ((generation + 1) * len(time_averages[param1])) == ((generation + 1) * len(time_averages[param2])) == len(pooled_ensemble[param1]) == len(pooled_ensemble[param2]):
                            print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
                            raise IOError('Sizes are not the same in a big way')
                        
                        # Get the variances of each type of series
                        ta_cov = ((generation + 1) * (
                                (time_averages[param1].copy() - time_averages.mean()[param1].copy()) * (time_averages[param2].copy() - time_averages.mean()[param2].copy()))).sum() / (
                                         len(pooled_ensemble) - 1)
                        pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
                        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                        
                        # Add it to the dataframe to output
                        out_temp[ind_out] = {
                            'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
                        }
                        
                        ind_out += 1
                        
                        if bootstraps != False:
                            for _ in np.arange(bootstraps):
                                # sample randomly from the distribution with replacement
                                time_averages_replacement = time_averages.sample(replace=True, frac=1)
                                
                                # Get the variances of each type of series
                                ta_cov = ((generation + 1) * (
                                        (time_averages_replacement[param1].copy() - time_averages_replacement.mean()[param1].copy()) * (
                                        time_averages_replacement[param2].copy() - time_averages_replacement.mean()[param2].copy()))).sum() / (
                                                 len(pooled_ensemble) - 1)
                                pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
                                
                                # Calculate the variance of the time-averages
                                gamma_ta_cov = ta_cov / pool_var
                                
                                # Add it
                                bs_temp[ind_bs] = {
                                    'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
                                }
                                
                                ind_bs += 1
                
                # To not repeat the calculation twice
                repeat.append(param1)
        
        # Append them to the big dataframes
        per_gen_ebp = per_gen_ebp.append(pd.DataFrame.from_dict(out_temp, "index"), ignore_index=True)
        per_gen_ebp_bs = per_gen_ebp_bs.append(pd.DataFrame.from_dict(bs_temp, "index"), ignore_index=True)
    
    # Important
    if per_gen_ebp.isnull().values.any():
        print('ERROR!')
        print(per_gen_ebp)
        exit()
    if per_gen_ebp_bs.isnull().values.any():
        print('ERROR! BS')
        print(per_gen_ebp_bs)
        exit()
    
    return [per_gen_ebp, per_gen_ebp_bs, total_length_df, total_length_bs_df]


def get_total_length_ebp(df, phenotypic_variables, expanding_mean, label, total_length_df, total_length_bs_df, bootstraps=False):
    # create the info dataframe with values being the cumulative mean of the lineage
    for lin_id in df.lineage_ID.unique():
        
        # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
        lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
        
        # add its time-average up until and including the total number of generations without NaNs
        to_add = pd.DataFrame.from_dict({
            'label': [label for _ in range(len(lineage))],
            'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
            'generation': [generation for generation in np.arange(len(lineage), dtype=int)]
        }).reset_index(drop=True)
        
        # Get the expanding mean and the cumulative sum, make sure they are robust to the first NaN in division ratio
        to_add_expanding_mean = pd.concat([lineage[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
        expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
        
        if to_add_expanding_mean.isnull().values.any():
            raise IOError('Got NaNs in the expanding mean, meaning something is wrong')
        
        lineage['generation'] = np.arange(len(lineage), dtype=int)
        
        if len(to_add_expanding_mean) != len(lineage):
            
            print(lineage)
            print(to_add_expanding_mean)
            raise IOError('Wronfg lenghts')
        
        if not np.array_equal(to_add_expanding_mean['generation'].values, lineage['generation'].values):
            print(lin_id)
            print(type(lineage['generation'].iloc[0]))
            print(type(to_add_expanding_mean['generation'].iloc[0]))
            raise IOError('Wronfg values of generation')
    
    # This is important
    assert not expanding_mean.isnull().values.any()
    
    # Get the EBP for the total length of the lineage
    time_averages = get_time_averages_df(df, phenotypic_variables)
    
    # for each parameter combination decompose the pooled variance
    tl_temp = {}
    tl_bs_temp = {}
    ind_tl_bs, ind_tl = 0, 0
    repeat = []
    for param1 in phenotypic_variables:
        for param2 in phenotypic_variables:
            if param2 not in repeat:
                # Get the variances of each type of series
                ta_cov = np.cov(time_averages[param1], time_averages[param2])[0, 1]
                
                pool_var = df[param1].std() * df[param2].std()
                
                # Calculate the variance of the time-averages
                gamma_ta_cov = ta_cov / pool_var
                
                # Add it to the dataframe to output
                tl_temp[ind_tl] = {
                    'param1': param1, 'param2': param2, 'n': len(time_averages[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label
                }
                
                # Continue the index of the dictionary we will make into a dataframe and add to what we return
                ind_tl += 1
                
                # If we decide to bootstrap here is where we do it.
                if bootstraps != False:
                    
                    for _ in np.arange(bootstraps):
                        # sample randomly from the distribution with replacement
                        time_averages_replacement = time_averages.sample(replace=True, frac=1).copy()
                        
                        # Get the covariance of each type of series
                        ta_cov = np.cov(time_averages[param1], time_averages[param2])[0, 1]
                        
                        pool_var = df[param1].std() * df[param2].std()
                        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                        
                        # Add it
                        tl_bs_temp[ind_tl_bs] = {
                            'param1': param1, 'param2': param2, 'n': len(time_averages[param1].unique()), 'generation': 'NA', 'gamma_ta': gamma_ta_cov, 'label': label
                        }
                        
                        # Continue the index of the dictionary we will make into a dataframe and add to what we return
                        ind_tl_bs += 1
        
        # To not repeat the calculation twice
        repeat.append(param1)
        
        # Append them to the big dataframes
        total_length_df = total_length_df.append(pd.DataFrame.from_dict(tl_temp, "index"), ignore_index=True)
        total_length_bs_df = total_length_bs_df.append(pd.DataFrame.from_dict(tl_bs_temp, "index"), ignore_index=True)
    
    # Important
    if total_length_df.isnull().values.any():
        print('ERROR! total_length_df')
        print(total_length_df)
        exit()
    if total_length_bs_df.isnull().values.any():
        print('ERROR! total_length_bs_df')
        print(total_length_bs_df)
        exit()
    
    return [total_length_df, total_length_bs_df]


def get_ebp_per_gen(df, phenotypic_variables, expanding_mean, label, per_gen_ebp, per_gen_ebp_bs, first_generation, final_generation, bootstraps=False):
    # For the pooled ensemble
    pooled_ensemble_total = pd.DataFrame(columns=df.columns)
    
    # create the info dataframe with values being the cumulative mean of the lineage
    for lin_id in df.lineage_ID.unique():
        
        # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
        lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
        
        # add its time-average up until and including the total number of generations without NaNs
        to_add = pd.DataFrame.from_dict({
            'label': [label for _ in range(len(lineage))],
            'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
            'generation': [generation for generation in np.arange(len(lineage), dtype=int)]
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
            raise IOError('Wronfg values of generation')
    
    # This is important
    assert not expanding_mean.isnull().values.any()
    
    # This is important
    assert not pooled_ensemble_total.isnull().values.any()
    
    # Calculate the ergodicity breaking parameter over all lineages in the dataset per generation
    for generation in np.arange(first_generation, final_generation + 1):
        print(generation)
        
        out_temp = {}
        bs_temp = {}
        ind_bs, ind_out = 0, 0
        
        # Get the time-averages of the lineages that have the generation we are looking for
        time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)].copy()  # Excludes all lineages that do not reach this generation
        
        # Get the array of which lineages are long enough, ie. have the amount of cycles we need
        long_enough_lineages = time_averages['lineage_ID'].unique()
        
        # Define the pooled ensemble of all the lineages that have at least this generation
        pooled_ensemble = pooled_ensemble_total[(pooled_ensemble_total['generation'] <= generation) & (pooled_ensemble_total['lineage_ID'].isin(long_enough_lineages))].copy()
        
        # That way we have a robust variance across lineages
        if len(long_enough_lineages) > 1:
            
            # for each parameter combination decompose the pooled variance
            repeat = []
            for param1 in phenotypic_variables:
                for param2 in phenotypic_variables:
                    if param2 not in repeat:
                        
                        # The mean is a linear operator
                        if np.abs(pooled_ensemble[param2].mean() - time_averages[param2].mean()) > 0.000001:
                            print(pooled_ensemble[param2].mean(), time_averages[param2].mean())
                            print(pooled_ensemble[param2].mean() - time_averages[param2].mean())
                            raise IOError('Means are not the same in a big way')
                        
                        # We must have the same amount of
                        if not ((generation + 1) * len(time_averages[param1])) == ((generation + 1) * len(time_averages[param2])) == len(pooled_ensemble[param1]) == len(pooled_ensemble[param2]):
                            print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
                            raise IOError('Sizes are not the same in a big way')
                        
                        # Get the variances of each type of series
                        ta_cov = ((generation + 1) * (
                                (time_averages[param1].copy() - time_averages.mean()[param1].copy()) * (time_averages[param2].copy() - time_averages.mean()[param2].copy()))).sum() / (
                                         len(pooled_ensemble) - 1)
                        pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
                        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                        
                        # Add it to the dataframe to output
                        out_temp[ind_out] = {
                            'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
                        }
                        
                        ind_out += 1
                        
                        if bootstraps != False:
                            for _ in np.arange(bootstraps):
                                # sample randomly from the distribution with replacement
                                time_averages_replacement = time_averages.sample(replace=True, frac=1)
                                
                                # Get the variances of each type of series
                                ta_cov = ((generation + 1) * (
                                        (time_averages_replacement[param1].copy() - time_averages_replacement.mean()[param1].copy()) * (
                                        time_averages_replacement[param2].copy() - time_averages_replacement.mean()[param2].copy()))).sum() / (
                                                 len(pooled_ensemble) - 1)
                                pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
                                
                                # Calculate the variance of the time-averages
                                gamma_ta_cov = ta_cov / pool_var
                                
                                # Add it
                                bs_temp[ind_bs] = {
                                    'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
                                }
                                
                                ind_bs += 1
                
                # To not repeat the calculation twice
                repeat.append(param1)
        
        # Append them to the big dataframes
        per_gen_ebp = per_gen_ebp.append(pd.DataFrame.from_dict(out_temp, "index"), ignore_index=True)
        per_gen_ebp_bs = per_gen_ebp_bs.append(pd.DataFrame.from_dict(bs_temp, "index"), ignore_index=True)
    
    # Important
    if per_gen_ebp.isnull().values.any():
        print('ERROR!')
        print(per_gen_ebp)
        exit()
    if per_gen_ebp_bs.isnull().values.any():
        print('ERROR! BS')
        print(per_gen_ebp_bs)
        exit()
    
    return [per_gen_ebp, per_gen_ebp_bs]


""" plot the ergodicity breaking parameter per generation (can be loglog or not) """



def plot_per_gen_ebp(df_dict, folder_name):
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
                    plt.ylabel(r'$\Gamma_n(${}$)$'.format(symbols['time_averages'][param1]))
                else:
                    plt.ylabel(r'$\Gamma_n(${}$, ${}$)$'.format(symbols['time_averages'][param1], symbols['time_averages'][param2]))
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


def plot_total_length_ebp(eb_df, symbolss, title):
    # The latex labels instead of the variable names
    eb_df = eb_df.replace(symbolss)
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    _, ax = plt.subplots(tight_layout=True, figsize=[6, 6])
    
    sns.barplot(x='param1', y='gamma_ta', data=eb_df[eb_df['param1'] == eb_df['param2']], hue='label', order=list(symbolss.values()), ax=ax,
                edgecolor='black')  # , label=r'$\overline{cov}$ \ $\sigma^2$'
    ax.yaxis.grid(True)
    # plt.ylim([0, .45])
    plt.title(title)
    ax.set_xlabel('')
    ax.set_ylabel(r'$\Gamma$')
    plt.legend(title='')
    # ax.get_legend(title='')
    plt.savefig('{}/EBP.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()


def main(args, first_generation, final_generation):
    # import/create the trace lineages
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    
    population_sampled = shuffle_info(physical_units, mm=args.MM)
    
    lineage_shuffled = shuffle_lineage_generations(physical_units, args.MM)
    
    # We keep the trap means here
    expanding_mean = pd.DataFrame(columns=['label', 'lineage_ID', 'generation'] + phenotypic_variables)
    
    # Where we keep the gammas
    per_gen_ebp = pd.DataFrame(columns=['param1', 'param2', 'n', 'generation', 'gamma_ta', 'label'])
    per_gen_ebp_bs = per_gen_ebp.copy()
    total_length_ebp = per_gen_ebp.copy()
    total_length_ebp_bs = per_gen_ebp.copy()
    
    # Calculate the cv and TA per lineage length
    for kind, df in zip(['Trace', 'Artificial', 'Shuffled'], [physical_units, population_sampled, lineage_shuffled]):
        
        print(kind)
        
        # Get the ebp for all lineages using the lengths we get from the measurements 
        total_length_ebp, total_length_ebp_bs = get_total_length_ebp(df, phenotypic_variables, expanding_mean, kind, total_length_ebp, total_length_ebp_bs, bootstraps=False)
        
        # Get the ebp for all lineages using the first n amount of generations, ie. an expanding window.
        # If there are other lineages that fall short we don't take them into account for that specific generation.
        # This means that the sample size for the number of lineages over which we compute the variance of the dataset can change between generations.
        # per_gen_ebp, per_gen_ebp_bs = get_ebp_per_gen(df, phenotypic_variables, expanding_mean, kind, per_gen_ebp, per_gen_ebp_bs, first_generation, final_generation, bootstraps=False)
        
        # per_gen_ebp, per_gen_ebp_bs, total_length_ebp, total_length_ebp_bs = get_both_per_gen_and_total_ebp(df, phenotypic_variables, expanding_mean, kind, per_gen_ebp, per_gen_ebp_bs, total_length_ebp, total_length_ebp_bs,
        #                                                                                         first_generation, final_generation, bootstraps=20)
        
        print(per_gen_ebp, per_gen_ebp_bs, total_length_ebp, total_length_ebp_bs, sep='\n' * 2)
        
        print('&' * 300)
    
    total_length_ebp.to_csv('gamma_ta_corrs_{}.csv'.format(args.data_origin), index=False)
    total_length_ebp_bs.to_csv('gamma_ta_corrs_bs_{}.csv'.format(args.data_origin), index=False)
    # per_gen_ebp.to_csv('gamma_ta_corrs_per_gen_{}.csv'.format(args.data_origin), index=False)
    # per_gen_ebp_bs.to_csv('gamma_ta_corrs_per_gen_bs_{}.csv'.format(args.data_origin), index=False)


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    total_ebp_dict = {}
    
    print(dataset_names)
    
    # Do all the Mother Machine data
    for data_origin in ['MC4100_25C', 'MC4100_27C', 'MC4100_37C']: # dataset_names:
        
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
        
        # In case they don't yet exist
        folder_name = args.figs_location + '/plot_per_gen_ebp'
        create_folder(args.figs_location)
        create_folder(folder_name)
        
        # Run the analysis
        main(args, first_generation=0, final_generation=30)
        
        # Import the analysis
        total = pd.read_csv('gamma_ta_corrs_{}.csv'.format(args.data_origin))
        
        # plot the total_length_ebp analysis for Trace and Shuffled
        plot_total_length_ebp(total[total['label'] != 'Shuffled'], symbols['time_averages'], title=data_origin)
        
        # # set a style on seaborn module
        # sns.set_context('paper')
        # sns.set_style("ticks", {'axes.grid': True})
        # _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
        # plot_total_length_ebp(, ax, symbols['time_averages'])
        # # save the figure
        # plt.savefig('{}/EBP.png'.format(args.figs_location), dpi=300)
        # # plt.show()
        # plt.close()
        
        print('*' * 200)
        exit()
