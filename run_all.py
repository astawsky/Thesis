
from DataframeCreation import mm_and_sm_processing, figure1_analysis, figure2_analysis, figure3_analysis, intergenerational_correlations, pair_gen_specific
from Plotting import figure1, figure2, figure3, intergenerational_figures, single_cell_correlations, growth_mechanism, histogram_of_lineage_lengths
from CustomFuncsAndVars.global_variables import mm_data_names, create_folder


def main():
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in mm_data_names + ['SM']:
        
        # Create the arguments for this function
        parser = argparse.ArgumentParser(description='Process Mother Machine and Sister Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
                            required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/RawData/' + data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                            required=False, default='population_lineages.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                            required=False, default='ergodicity_breaking_parameter.csv')
        parser.add_argument('-kld', '--kld', metavar='', type=str,
                            help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                            required=False, default='kullback_leibler_divergences.csv')
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Figures/' + data_origin)
        parser.add_argument('-scc', '--scc', metavar='', type=str, help='Where the single cell correlation figures are saved.',
                            required=False, default='single_cell_correlations')
        parser.add_argument('-scch', '--scch', metavar='', type=str, help='Where the single cell correlation heatmap figures are saved.',
                            required=False, default='single_cell_correlations_heatmaps')
        parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                            required=False, default='trace_centered.csv')
        
        # We need this for the SM data processing
        if data_origin == 'SM':
            parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
            parser.add_argument('-SL', '--sl_infiles', metavar='', type=str,
                                help='Location of Sister Lineage Raw Data', required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/RawData/SM/SL')
            parser.add_argument('-NL', '--nl_infiles', metavar='', type=str,
                                help='Location of Neighboring Lineage Raw Data', required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/RawData/SM/NL')
        else:
            parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        
        # Finalize the arguments
        args = parser.parse_args()
        
        # Make sure the folders where we place the data are created already
        create_folder(args.raw_data)
        create_folder(args.save_folder)
        
        if data_origin == 'SM':  # Get SM data
            # How long did it take to do the mother machine?
            mm_time = time.time() - first_time

            # Get SM data
            mm_and_sm_processing.main_sm(args)
        else:  # Get MM data
            mm_and_sm_processing.main_mm(args)

        # exit()

        # Calculate the dataframes necessary
        figure1_analysis.main(args)
        # Plot
        histogram_of_lineage_lengths.main(args)
        figure1.individuals(args)
        single_cell_correlations.main(args)
        growth_mechanism.main(args)
        
        print('*' * 200)
        # exit()
    
    # How much time does the SM data take?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))


if __name__ == '__main__':
    main()
