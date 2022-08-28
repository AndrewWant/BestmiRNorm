# Code for identifying normaliser combinations
# miRNA_normalisation
# Identification of an optimal normaliser combination

import numpy as np
import pandas as pd
import os
from tkinter import filedialog
import tkinter
from itertools import combinations

def get_filename(type="raw"):
    filename_args = {"raw": {"filetypes":(("Excel file","*.xlsx"),
                                          ("Excel file", "*.xls"))}
                     "normaliser": {"filetypes": "*.txt"}}
    cwd = os.getcwd()
    root = tkinter.Tk()
    filename = filedialog.askopenfilename(parent=root,
                                          initialdir=cwd,
                                          title="Please select your {} file".format(type),
                                          **filename_args[type])
    filename = filename.replace("/", "\\")
    root.destroy()
    return filename


def get_candidate_normalisers(normaliser_filename):
    with open(normaliser_filename, "r") as open_normalisers:
        normalisers = open_normalisers.readlines()
        normalisers = [x.strip("\n") for x in normalisers]
        print("You have selected {} normalisers".format(len(normalisers)))
        return normalisers

        
def read_xl(xl_filename, sheet_name):
    df = pd.read_excel(xl_filename,
                       sheet_name=sheet_name,
                       header=[0, 1, 2, 3],
                       index_col=[0])
    return df
    

def generate_normalisers(normaliser_locations=(0, 1, 3, 4, 5, 6)):
    """
    Generate combinations normaliser indices based on normaliser_locations, included optional kw for endogenous control
    location(s)
    :param normaliser_locations:
    :return:
    """
    for i in range(1, len(normaliser_locations) + 1):
        c = combinations(normaliser_locations, i)
        for x in c:
            yield np.asarray(x)


def log_fold_change_from_cq(cq_dataframe,
                            normalisers,
                            control="HC",
                            rtn_log=True):
    """
    A single function for generating the log2(fold change) from averaged Cq values
    :param cq_dataframe:
    :param normalisers:
    :param control:
    :param rtn_log:
    :return:
    """
    if len(normalisers) == 1:
        norm_mean = cq_dataframe.loc[normalisers, :].values
    else:
        norm_mean = cq_dataframe.loc[normalisers, :].mean()
    dcq_df = cq_dataframe - norm_mean
    ddcq_df = dcq_df - dcq_df.loc[:, (control, slice(None))].mean(axis=1).values.reshape((-1, 1))
    fold_change = 2 ** (ddcq_df * -1)
    if rtn_log:
        return np.log2(fold_change)
    else:
        return fold_change


def make_result_dataframe(cq_dataframe,
                          qty_normalisers
                          ):
    """
    Creation of an empty DataFrame that will enable more effective storage of experimental data.
    :param cq_dataframe:
    :param qty_normalisers: This is the number of combinations that will be generated!
    :return:
    """
    r_idx = pd.MultiIndex.from_product((range(qty_normalisers), cq_dataframe.index),
                                       names=("combination_code", "miRNA_name")
                                       )
    dataframe = pd.DataFrame(index=r_idx,
                             columns=cq_dataframe.columns,
                             dtype=float
                             )
    dataframe[:] = 0
    return dataframe


def populate_dataframe(cq_dataframe: pd.DataFrame,
                       output_dataframe: pd.DataFrame,
                       normaliser_generator):
    """
    Fills a target dataframe with log(fold change) or delta Cq data from a dictionary.
    :param cq_dataframe:
    :param output_dataframe:
    :param normaliser_generator:
    :return:
    """
    retained_normalisers = {}
    mirnas = cq_dataframe.index
    for index_, normalisers in enumerate(normaliser_generator):
        retained_normalisers[index_] = mirnas[normalisers]
        output_dataframe.loc[(index_, mirnas), :] = log_fold_change_from_cq(cq_dataframe,
                                                                            mirnas[normalisers]).values
    return output_dataframe, retained_normalisers


def analyse_logfc_data(log_fc_dataframe,
                       normalisers,
                      positive_class,
                      control_class):
    """
    Comparing log(FC) distributions based on population metrics, and creating combination metrics for AD and HC
    populations
    :param log_fc_dataframe:
    :param normalisers:
    :return summary_data:
    """
    ks_out = {}
    for idx in log_fc_dataframe.index.get_level_values(0):
        pos_data = log_fc_dataframe.loc[(idx, normalisers), (positive_class, slice(None))].values.flatten()
        cont_data = log_fc_dataframe.loc[(idx, normalisers), (control_class, slice(None))].values.flatten()
        k, p = ks_2samp(pos_data, cont_data)
        pos_mean = pos_data.mean()
        pos_std = pos_data.std()
        cont_mean = cont_data.mean()
        cont_std = cont_data.std()

        ks_out[idx] = (k, p, pos_mean, pos_std, cont_mean, cont_std)

    header = ("KS-score", 
              "p-value", 
              "{}_mean".format(positive_class), 
              "{}_std".format(positive_class), 
              "{}_mean".format(control_class), 
              "{}_std".format(control_class))
    # new_idx = pd.MultiIndex.from_arrays((range(len(ks_out)), ks_out.keys()))
    df = pd.DataFrame(ks_out.values(), index=ks_out.keys(), columns=header)
    df["root_square_sum_mean"] = np.sqrt(df[header[2]] ** 2 + df[header[4]] ** 2)
    df["mean_stds"] = df[[header[3], header[-1]]].mean(axis=1)

    summary_data = df.loc[:, np.asarray(("KS-score", "root_square_sum_mean", "mean_stds"))]

    return summary_data


def rank_weight_generator(weights=(1, 2, 3)):
    """Create spread of all rank weights for 1-3 times weighting

       Note: (2, 2, 2) and (3, 3, 3) are discarded, because they will be equivalent to (1, 1, 1)
    """
    for rank_weights in product(weights, repeat=3):
        i, j, k = rank_weights
        if (i > 1) and (i == j == k):
            continue
        else:
            yield rank_weights


def rank_data(stat_summary_df,
              rank_weights=(1, 1, 1),
              rank_basis="sum"):
    """
    Create a DataFrame of the rankings for each metric, with the sum of all metrics, including application of weights to
    metrics
    :param stat_summary_df:
    :param rank_weights:
    :param rank_basis:
    :return ranking_df:
    """
    indices = {}

    for col in ("KS-score", "root_square_sum_mean", "mean_stds"):
        sorted_data = stat_summary_df.sort_values(by=col, ascending=False)
        for rank, norm_code in enumerate(sorted_data.index):
            if norm_code in indices:
                indices[norm_code].append(rank)
            else:
                indices[norm_code] = [rank]

    ranking_df = pd.DataFrame(indices.values(),
                              index=indices.keys(),
                              columns=("KS-score", "root_square_sum_mean", "mean_stds"))
    ranking_df *= rank_weights
    ranking_df["sum"] = ranking_df.sum(axis=1)
    return ranking_df.sort_values(by=rank_basis, ascending=False)

                     
def generate_suppl_ranked(cq_dataframe,
                          weights,
                          normaliser_indices,
                          normalisers,
                          rank_basis="sum",
                          control="HC"):
    normaliser_conversion_dict = dict([(idx, np.asarray(cq_dataframe.index[normalisers])) for (idx, normalisers) in
                                       enumerate(generate_normalisers(normaliser_indices))])
    norm_generator = generate_normalisers(normaliser_indices)
    result_df = make_result_dataframe(cq_dataframe,
                                      len(normaliser_conversion_dict))
    combination_df, combination_normalisers = populate_dataframe(cq_dataframe,
                                                                 result_df,
                                                                 norm_generator)
    population_analysis = analyse_logfc_data(combination_df,
                                             normalisers,
                                            control=control)
    ranked_df = rank_data(population_analysis,
                          weights,
                          rank_basis)
    return ranked_df, normaliser_conversion_dict, combination_df, population_analysis


def compare_ranked_normalisers(ranked_df,
                               conversion_dict: dict,
                               list_of_normalisers,
                               top_n_normalisers=10):
    normaliser_number_value = dict([(ctrl, i) for (i, ctrl) in enumerate(list_of_normalisers, 1)])
    normaliser_number_value[""] = 0
    number_switched_norms = {}
    top_n = ranked_df.index[:top_n_normalisers].to_list()
    top_n_output = []
    for entry in top_n:
        normalisers = conversion_dict[entry]
        normaliser_with_blanks = []
        for ctrl in list_of_normalisers:
            if ctrl in normalisers:
                normaliser_with_blanks.append(normaliser_number_value[ctrl])
            else:
                normaliser_with_blanks.append(normaliser_number_value[""])
        top_n_output.append(normaliser_with_blanks)
    number_switched_norms[i] = top_n_output
    return number_switched_norms
       
    
