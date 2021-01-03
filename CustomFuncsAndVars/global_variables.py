#!/usr/bin/env bash


""" create a folder with the name of the string inputted """


def create_folder(filename):
    import os
    
    # create the folder if not created already
    try:
        # Create target Directory
        os.mkdir(filename)
        print("Directory", filename, "Created ")
    except FileExistsError:
        pass


""" Returns dataframe with same size containing the time-averages of each phenotypic variable instead of the local value """


def get_time_averages_df(info, phenotypic_variables):  # MM
    import numpy as np
    
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


""" shuffle the generations inside a lineage to destroy inter-generational correlations """


def shuffle_lineage_generations(df, mm):
    import pandas as pd
    import numpy as np
    
    if mm:
        
        # This is where we will keep all the data and later save as a csv
        no_epi_df = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
        
        for lin_id in df.lineage_ID.unique():
            lineage = df[(df['lineage_ID'] == lin_id)]
            
            new_lineage = lineage[phenotypic_variables].sample(frac=1, replace=False)  # .reset_index(drop=True)
            new_lineage['lineage_ID'] = lineage.lineage_ID.unique()[0]
            new_lineage[phenotypic_variables] = lineage[phenotypic_variables].sample(frac=1, replace=False)
            new_lineage['generation'] = np.arange(len(lineage))
            
            # adds it to the dataframe
            no_epi_df = no_epi_df.append(new_lineage, ignore_index=True)
    else:
        
        # This is where we will keep all the data and later save as a csv
        no_epi_df = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'dataset', 'trap_ID', 'trace', 'generation'])
        
        # loop for every dataset and its traps and traces
        for dataset in np.unique(df['dataset']):
            
            for trap_id in np.unique(df[(df['dataset'] == dataset)]['trap_ID']):
                
                for trace in ['A', 'B']:
                    
                    lineage = df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)]
                    
                    new_lineage = lineage[phenotypic_variables].sample(frac=1, replace=False)  # .reset_index(drop=True)
                    new_lineage['lineage_ID'] = lineage.lineage_ID.unique()[0]
                    new_lineage['dataset'] = dataset
                    new_lineage['trap_ID'] = trap_id
                    new_lineage['trace'] = trace
                    new_lineage['generation'] = np.arange(len(lineage))
                    
                    # adds it to the dataframe
                    no_epi_df = no_epi_df.append(new_lineage, ignore_index=True)
                
                # exit()
                #
                # # With these we mask the dataframes to get the values we want, mainly a lineage's empirical time-average and standard deviation
                # trace_mask = (df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)
                #
                # # how many generations we have in this lineage
                # max_gen = len(df[trace_mask])
                #
                # # What we are going to add to the
                # to_add = {
                #     'dataset': [dataset for _ in range(max_gen)],
                #     'trap_ID': [trap_id for _ in range(max_gen)],
                #     'trace': [trace for _ in range(max_gen)],
                #     'generation': np.arange(max_gen)
                # }
                #
                # # shuffle the lineage generations so epigenetic MD memory is lessened and or destroyed
                # shuffled = df[trace_mask].sample(frac=1, replace=False).reset_index(drop=True)
                #
                # # creates the iid array of samples from every cycle variableeter lineage distribution
                # for variable in phenotypic_variables:
                #     to_add.update({variable: shuffled[variable]})
                #
                # # adds it to the dataframe
                # no_epi_df = no_epi_df.append(pd.DataFrame(to_add), ignore_index=True)
    
    # Reset the index as a matter of cleanliness since we are not gonna work with it anyways
    return no_epi_df.reset_index(drop=True)


""" truncate all lineages in a dataframe to the min_gens parameter """


def limit_lineage_length(df, min_gens):
    import numpy as np
    
    new_df = pd.DataFrame(columns=df.columns)
    for dataset in df['dataset'].unique():
        for trap_id in df[(df['dataset'] == dataset)]['trap_ID'].unique():
            for trace in ['A', 'B']:
                lineage = df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)].copy()
                
                if len(lineage) >= min_gens:  # It has at least these amount of generations and we will truncate them there
                    # (min_gens - 1) is because the first generation is equal to zero, not 1
                    new_lineage = lineage[lineage['generation'] <= (min_gens - 1)][phenotypic_variables].sample(frac=1, replace=False)  # .reset_index(drop=True)
                    new_lineage['dataset'] = dataset
                    new_lineage['trap_ID'] = trap_id
                    new_lineage['trace'] = trace
                    new_lineage['generation'] = np.arange(min_gens)
                    
                    # adds it to the dataframe
                    new_df = new_df.append(new_lineage, ignore_index=True)
                else:
                    pass
    
    return new_df


""" Feed in a dataframe of our format and subtract the time-averages of each lineage """


def trace_center_a_dataframe(df):
    import pandas as pd
    import numpy as np
    
    tc_df = pd.DataFrame(columns=df.columns)
    
    for dataset in df['dataset'].unique():
        for trap_id in df[df['dataset'] == dataset]['trap_ID'].unique():
            for trace in ['A', 'B']:
                lineage = df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)]
                
                new_lineage = lineage[phenotypic_variables] - lineage[phenotypic_variables].mean()
                new_lineage['dataset'] = dataset
                new_lineage['trap_ID'] = trap_id
                new_lineage['trace'] = trace
                new_lineage['generation'] = np.arange(len(lineage))
                
                tc_df = tc_df.append(new_lineage, ignore_index=True)
    
    return tc_df


""" Create all information dataframe where the lineage lengths are kept constant but the cells in the trace itself are randomly sampled from the population without replacement """


def shuffle_info_OLD(info):
    import pandas as pd
    import numpy as np
    
    # Give it a name, contains S, NL
    new_info = pd.DataFrame(columns=info.columns)
    
    # what is the trace length of each trace? This is the only thing that stays the same
    sizes = {'{} {} {}'.format(dataset, trap_ID, trace): len(info[(info['trap_ID'] == trap_ID) & (info['trace'] == trace) & (info['dataset'] == dataset)]) for dataset in np.unique(info['dataset']) for
             trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']) for trace in ['A', 'B']}
    
    for dataset in np.unique(info['dataset']):
        for trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']):
            for trace in ['A', 'B']:
                # trace length
                size = sizes['{} {} {}'.format(dataset, trap_ID, trace)]
                
                # sample from the old dataframe
                samples = info.sample(replace=False, n=size)
                
                # drop what we sampled
                info = info.drop(index=samples.index)
                
                # add some correct labels even though they don't matter so that the dataframe structure is still intact
                samples['trap_ID'] = trap_ID
                samples['trace'] = trace
                samples['generation'] = np.arange(size)
                samples['dataset'] = dataset
                
                # add them to the new, shuffled dataframe
                new_info = new_info.append(samples, ignore_index=True)
    
    return new_info


""" Create all information dataframe where the lineage lengths are kept constant but the cells in the trace itself are randomly sampled from the population without replacement """


def shuffle_info(info, mm):
    import pandas as pd
    import numpy as np
    
    # Give it a name, contains S, NL
    new_info = pd.DataFrame(columns=info.columns)
    
    # what is the trace length of each trace? This is the only thing that stays the same
    if mm:
        sizes = {'{}'.format(lineage_ID): len(info[(info['lineage_ID'] == lineage_ID)]) for lineage_ID in info['lineage_ID'].unique()}
        
        for lineage_id in np.unique(info['lineage_ID']):
            # trace length
            size = sizes['{}'.format(lineage_id)]
            
            # sample from the old dataframe
            samples = info.sample(replace=False, n=size)
            
            # drop what we sampled
            info = info.drop(index=samples.index)
            
            # add some correct labels even though they don't matter so that the dataframe structure is still intact
            samples['lineage_ID'] = lineage_id
            samples['generation'] = np.arange(size)
            
            # add them to the new, shuffled dataframe
            new_info = new_info.append(samples, ignore_index=True)
    else:
        sizes = {'{} {} {}'.format(dataset, trap_ID, trace): len(info[(info['trap_ID'] == trap_ID) & (info['trace'] == trace) & (info['dataset'] == dataset)]) for dataset in np.unique(info['dataset'])
                 for
                 trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']) for trace in ['A', 'B']}
        
        lineage_id = 0
        for dataset in np.unique(info['dataset']):
            for trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']):
                for trace in ['A', 'B']:
                    # trace length
                    size = sizes['{} {} {}'.format(dataset, trap_ID, trace)]
                    
                    # sample from the old dataframe
                    samples = info.sample(replace=False, n=size)
                    
                    # drop what we sampled
                    info = info.drop(index=samples.index)
                    
                    # add some correct labels even though they don't matter so that the dataframe structure is still intact
                    samples['dataset'] = dataset
                    samples['trap_ID'] = trap_ID
                    samples['trace'] = trace
                    samples['lineage_ID'] = lineage_id
                    samples['generation'] = np.arange(size)
                    
                    # add them to the new, shuffled dataframe
                    new_info = new_info.append(samples, ignore_index=True)
                    
                    # Go on to the next ID
                    lineage_id += 1
    
    return new_info


""" What we use for all main-text figures """


def seaborn_preamble():
    import seaborn as sns
    
    # stylistic reasons
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})


""" set of bootstrapped pearson correlations """


def boot_pearson(x, y, number_of_straps):
    import numpy as np
    from random import choices
    from scipy.stats import pearsonr
    x = np.array(x)
    y = np.array(y)
    r_array = []
    for n in range(number_of_straps):
        indices = choices(np.arange(len(x)), k=len(x))
        r = pearsonr(x[indices], y[indices])[0]
        r_array.append(r)
    
    return np.array(r_array)


""" We use this to cut any lineage-processed dataframe, ie. '.../phyiscal_units.csv' and '.../trace_centered.csv' so that the A and B pairs have the same lineage length """


def cut_uneven_pairs(info):
    import numpy as np
    import pandas as pd
    
    # Correct for traces that aren't the same size
    old_info = info.copy()
    new_info = pd.DataFrame(columns=info.columns)
    
    for dataset in np.unique(old_info['dataset']):
        for trap_id in np.unique(old_info[(old_info['dataset'] == dataset)]['trap_ID']):
            # print(trap_id)
            min_gen = min(max(old_info[(old_info['trap_ID'] == trap_id) & (old_info['trace'] == 'A') & (old_info['dataset'] == dataset)]['generation']),
                          max(old_info[(old_info['trap_ID'] == trap_id) & (old_info['trace'] == 'B') & (old_info['dataset'] == dataset)]['generation']))
            
            # What we will add
            a_to_add = old_info[(old_info['trap_ID'] == trap_id) & (old_info['trace'] == 'A') & (old_info['dataset'] == dataset) & (old_info['generation'] <= min_gen)]
            b_to_add = old_info[(old_info['trap_ID'] == trap_id) & (old_info['trace'] == 'B') & (old_info['dataset'] == dataset) & (old_info['generation'] <= min_gen)]
            
            # Check that they are the same size
            if len(a_to_add) != len(b_to_add):
                print(a_to_add)
                print(b_to_add)
                raise IOError()
            
            # save them to a new dataframe
            new_info = new_info.append(a_to_add, ignore_index=True)
            new_info = new_info.append(b_to_add, ignore_index=True)
        
        # good practice
        new_info = new_info.reset_index(drop=True)
    
    return new_info



import pandas as pd
import seaborn as sns

phenotypic_variables = ['div_then_fold', 'div_and_fold', 'fold_then_div', 'fold_growth', 'division_ratio', 'added_length', 'generationtime', 'length_birth', 'length_final', 'growth_rate']
symbols = {
    'physical_units': dict(
        zip(phenotypic_variables, [r'$f_n e^{\phi_{n+1}}$', r'$f_n e^{\phi_{n}}$', r'$f_{n+1} e^{\phi_n}$', r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'])),
    'trace_centered': dict(zip(phenotypic_variables,
                               [r'$\delta (f_n e^{\phi_{n+1}})$', r'$\delta (f_n e^{\phi_{n}})$', r'$\delta (f_{n+1} e^{\phi_n})$', r'$\delta \phi$', r'$\delta f$', r'$\delta \Delta$',
                                r'$\delta \tau$', r'$\delta x_0$', r'$\delta x_\tau$', r'$\delta \alpha$'])),
    'time_averages': dict(
        zip(phenotypic_variables,
            [r'$\overline{f_n e^{\phi_{n+1}}}$', r'$\overline{f_n e^{\phi_{n}}}$', r'$\overline{f_{n+1} e^{\phi_n}}$', r'$\overline{\phi}$', r'$\overline{f}$', r'$\overline{\Delta}$',
             r'$\overline{\tau}$', r'$\overline{x_0}$', r'$\overline{x_\tau}$', r'$\overline{\alpha}$'])),
    'model': dict(zip(phenotypic_variables, [r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'])),
    'old_var': dict(zip(phenotypic_variables, [r'$\phi_{old}$', r'$f_{old}$', r'$\Delta_{old}$', r'$\tau_{old}$', r'$x_{0, old}$', r'$x_{\tau, old}$', r'$\alpha_{old}$'])),
    'young_var': dict(zip(phenotypic_variables, [r'$\phi_{young}$', r'$f_{young}$', r'$\Delta_{young}$', r'$\tau_{young}$', r'$x_{0, young}$', r'$x_{\tau, young}$', r'$\alpha_{young}$'])),
    'with_n_ta': dict(zip(phenotypic_variables, [r'$\overline{\phi_{n}}$', r'$\overline{f_{n}}$', r'$\overline{\Delta_{n}}$', r'$\overline{\tau_{n}}$', r'$\overline{x_{0, n}}$',
                                                 r'$\overline{x_{\tau, n}}$', r'$\overline{\alpha_{n}}$'])),
    'unique_ta': dict(
        zip(phenotypic_variables,
            [r'$\overline{f_n e^{\phi_{n+1}}}$', r'$\overline{f_n e^{\phi_{n}}}$', r'$\overline{f_{n+1} e^{\phi_n}}$', r'$\overline{\phi}$', r'$\overline{f}$', r'$\overline{\Delta}$',
             r'$\overline{\tau}$', r'$\overline{x_0}$', r'$\overline{x_\tau}$', r'$\overline{\alpha}$'])),
    # 'new_pu': dict(zip(phenotypic_variables,
    #                    [r'$\overline{\phi}$', r'$\overline{f}$', r'$\overline{\Delta}$', r'$\overline{\tau}$', r'$\overline{x_0}$', r'$\overline{x_\tau}$', r'$\overline{\alpha}$',
    #                     r'$\overline{\frac{x_{0, n+1}}{x_{0, n}}}$']))
}
datasets = ['SL', 'NL', 'CTRL']
units = dict(zip(phenotypic_variables, [r'', r'', r'', r'', r'', r'$\mu m$', r'$(hr)$', r'$(\mu m)$', r'$(\mu m)$', r'$(\frac{1}{hr})$']))
hierarchy = [r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$']
# This is a manual thing to make the distribution look better without the outliers
bounds = {
    'generationtime': [0.01, 1], 'fold_growth': [.27, 1], 'growth_rate': [.8, 1.8], 'length_birth': [.7, 4.5], 'length_final': [2.3, 8.3], 'division_ratio': [.35, .65], 'added_length': [.4, 5]
}
symbols_bounds = {symbols['physical_units'][key]: val for key, val in bounds.items()}
dataset_names = ['1015_NL', '062718_SL', '071318_SL', '072818_SL_NL', '101218_SL_NL', 'Pooled_SM', '8-31-16 Continue', 'MC4100_25C', 'MC4100_27C', 'MC4100_37C']  # , 'LAC_M9', , '4-28-2017'
# dataset_names2 = ['SM', 'MG1655_inLB_LongTraces', 'Maryam_LongTraces', 'lambda_LB', 'MC4100_25C', 'MC4100_27C', 'MC4100_37C']  # , 'LAC_M9', , '4-28-2017'
cmap = sns.color_palette('tab10')
