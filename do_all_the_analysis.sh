#!/usr/bin/env bash


echo -----------------------------------

# Create the physical_units and trace_centered dataframes in .csv format in Data folder
python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/lineage_data_processing.py -h

python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/lineage_data_processing.py

echo -----------------------------------

# Create the physical_units and trace_centered dataframes in .csv format in Data folder with Control dataset
python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/pair_data_processing.py -h

python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/pair_data_processing.py

echo -----------------------------------



echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure1_analysis.py -h

python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure1_analysis.py

echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure1.py -h

python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure1.py

echo -----------------------------------



echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure2_analysis.py -h

python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure2_analysis.py

echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure2.py -h

python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure2.py

echo -----------------------------------



echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure3_analysis.py -h

python /Users/alestawsky/PycharmProjects/Thesis/DataframeCreation/figure3_analysis.py

echo -----------------------------------

# Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.
python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure3.py -h

python /Users/alestawsky/PycharmProjects/Thesis/Plotting/figure3.py

echo -----------------------------------
