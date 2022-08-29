import miRNA_normalisers as miR_norm
import check_input_files as chk_files

import numpy as np
import os
import pandas as pd
from datetime import datetime
import sys


def find_negative_class_label(cq_dataframe, positive_class_label):
    available_classes = np.unique(cq_dataframe.columns.get_level_values("Biological Group"))
    negative_class = [x for x in available_classes if x != positive_class_label]
    return negative_class[0]


def ask_rank_weights():
    rank_weights = []
    score_components = ("Kolmogorov-Smirnov score", "Mean distance from zero", "Median of standard deviations")
    for component in score_components:
        next_rank_weight = input("Please indicate weight to apply to {}: (numerical values only) ".format(component))
        try:
            weight_as_number = float(next_rank_weight)
            rank_weights.append(weight_as_number)
        except ValueError:
            rank_weights.append(next_rank_weight)
    if np.any(rank_weights == None):
        print("One or more weights were not valid, your weights were: {}, {}, {}\nApplying equal weight...".format(*rank_weights))
        return (1, 1, 1)
    else:
        return rank_weights
    

def create_output_filepath(input_filepath):
    folder_path = os.sep.join(input_filepath.split(os.sep)[:-1])
    return folder_path

def create_output_filename(label):
    now = datetime.now()
    now_string = now.strftime("%Y-%m-%d_%H%M%S")
    return "_".join((now_string, label, "output_file"))


xl_file = miR_norm.get_filename()
sheet_name = input("Enter Sheet Name with Your Data, or press enter to use the default from the template: ")
positive_class = input("Enter the label of your positive class (e.g. disease-bearing subjects): ")

# if sheet_name == "":
sheet_name == "Results"
normaliser_file = miR_norm.get_filename("normaliser")

normalisers = miR_norm.get_candidate_normalisers(normaliser_file)

if chk_files.all_checks(miR_norm.read_xl(xl_file, sheet_name), normaliser_file, positive_class):
    print("All file checks passed: Proceeding to analysis")
else:
    print("Correct files where indicated, and re-run")
    sys.exit()

specified_weights = ask_rank_weights()
miRNA_data = miR_norm.read_xl(xl_file, sheet_name)
negative_class = find_negative_class_label(miRNA_data, positive_class)
normaliser_locations = chk_files.get_normaliser_locations(miRNA_data, normalisers)
normaliser_gen = miR_norm.generate_normalisers(normaliser_locations)
number_of_normaliser_combinations = (2 ** len(normalisers)) - 1
ranked_df, normaliser_conversion_dict, combination_df, population_analysis = miR_norm.generate_suppl_ranked(miRNA_data, 
                                                                                                            specified_weights,
                                                                                                            normaliser_locations,
                                                                                                            normalisers,
                                                                                                            control=negative_class)


output_folder = create_output_filepath(xl_file)
output_dict = {"ranked_df": ranked_df,
               "normaliser_conversion_dict": normaliser_conversion_dict,
               "combination_df": combination_df,
               "population_analysis": population_analysis}
                     
for output_label, output in output_dict.items():
    output_filename = create_output_filename(output_label)
    output_filepath = os.sep.join((output_folder, output_filename))
    if isinstance(output, pd.DataFrame):
        output_filename = output_filepath + ".xlsx"
        output.to_excel(output_filename)
    elif isinstance(output, dict):
        output_filename = output_filepath + ".txt"
        with open(output_filename, "a") as opened_file:
            for key, value in output.items():
                value_output = "{}" * len(value)
                # value_output = value_output.format(value)
                opened_file.write("{}, {}, {}".format(key, "\t", value))
            
