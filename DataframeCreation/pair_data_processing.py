#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, cut_uneven_pairs
import pandas as pd
import numpy as np
import argparse
from random import sample, seed
from itertools import combinations


""" Add Control, cut extra generations, create the trace centered and time average dataframes and save them """


def add_control(info, number_of_control_lineage_pairs):
    
    """
    Here we will create the Control dataset and equal the length of pair lineages.
    """
    
    # For the Control
    traces_to_choose_from = info['dataset'] + ', ' + info['trap_ID'].astype(str)
    
    # Choose randomly combinations of trap IDs for the A and B traces
    possible_combos = sample(list(combinations(np.unique(traces_to_choose_from), 2)), number_of_control_lineage_pairs)
    
    # start them as new IDs
    new_id = max(np.unique(info['trap_ID'])) + 1
    
    # loop through all traces paired together for Control
    for combo in possible_combos:
        # define the trace ID
        dataset_a, id_a = combo[0].split(', ')
        dataset_b, id_b = combo[1].split(', ')
        
        # define the trace
        a_trace = info[(info['trap_ID'] == int(id_a)) & (info['trace'] == 'A') & (info['dataset'] == dataset_a)].copy()
        # change the trap id to a new id number even though it is a copy of the other ID number
        a_trace['trap_ID'] = new_id
        # Change the dataset to Control
        a_trace['dataset'] = 'CTRL'
        # Add it to the information dataframe
        info = pd.concat([info, a_trace], axis=0)
        
        # Do the same
        b_trace = info[(info['trap_ID'] == int(id_b)) & (info['trace'] == 'B') & (info['dataset'] == dataset_b)].copy()
        b_trace['trap_ID'] = new_id
        b_trace['dataset'] = 'CTRL'
        info = pd.concat([info, b_trace], axis=0)
        
        # create a new ID number to add in the next loop
        new_id += 1
    
    # Just for my preference even though we don't rely on indices
    info = info.reset_index(drop=True)
    
    return info


def trace_centered_and_trace_means():
    
    # copy it just in case
    info = physical_units.copy()
    
    """ Here we get the trace centered units and the trace time averages """
    
    # read the csv file where we keep the data
    trace_centered = info.copy()
    
    # We keep the trap means here
    trace_means_df = pd.DataFrame(columns=['dataset', 'trap_ID', 'trace', 'max_gen', 'generation'] + phenotypic_variables)
    
    # specify a lineage
    for dataset in np.unique(info['dataset']):
        for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
            for trace in ['A', 'B']:
                
                # the values of the lineage we get from physical units
                lineage = info[(info['trap_ID'] == trap_id) & (info['dataset'] == dataset) & (info['trace'] == trace)].copy()
                
                # Trace centered
                trace_centered.loc[lineage.index, phenotypic_variables] = lineage[phenotypic_variables] - lineage[phenotypic_variables].mean()
                
                # add its time-average
                to_add = {
                    'dataset': [dataset for _ in np.arange(len(lineage))],
                    'trap_ID': [trap_id for _ in np.arange(len(lineage))],
                    'trace': [trace for _ in np.arange(len(lineage))],
                    'max_gen': [len(lineage) for _ in np.arange(len(lineage))],
                    'generation': np.arange(len(lineage))
                }
                to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables})
                to_add = pd.DataFrame(to_add)
                trace_means_df = trace_means_df.append(to_add, ignore_index=True)

    trace_means_df = trace_means_df[phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'max_gen', 'generation']]
    
    assert len(trace_centered) == len(info) == len(trace_means_df)
    
    return [trace_centered, trace_means_df]
    
    


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
# parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default='physical_units_with_control.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default='time_averages_with_control.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default='trace_centered_with_control.csv')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-ctrl_number', '--number_of_control_lineage_pairs', metavar='', type=int, help='How many Control lineage pairs should we make.',
                    required=False, default=85)
parser.add_argument('-seed', '--random_seed', metavar='', type=int, help='Seed for python randomness generator.',
                    required=False, default=42)


args = parser.parse_args()


create_folder(args.save_folder)


# So we can compare the trace-centered and the physical units
seed(args.random_seed)

# read the csv file where we keep the data
from_lineage_processing = pd.read_csv('{}/physical_units.csv'.format(args.save_folder))

# add the control dataset to the measured and labeled data and cut equalize the lineage length betwen two pair lineages by using the least longest length between the lineages.
physical_units = add_control(from_lineage_processing, args.number_of_control_lineage_pairs)
physical_units = cut_uneven_pairs(physical_units)

# create the trace_centered and trace_means
trace_centered, trace_means_df = trace_centered_and_trace_means()

# save them to the Data folder
physical_units.to_csv('{}/{}'.format(args.save_folder, args.pu), index=False)
trace_means_df.to_csv('{}/{}'.format(args.save_folder, args.ta), index=False)
trace_centered.to_csv('{}/{}'.format(args.save_folder, args.tc), index=False)
