import pandas as pd
import numpy as np

def check_xl_column_names(xl_dataframe):
    canon_names = ("Biological Group", "subject", "sex", "age")
    if np.all(xl_dataframe.columns.names == canon_names):
        return True
    else:
        for idx, (xl_name, canon_name) in enumerate(zip(xl_dataframe.columns.names, canon_names)):
            if xl_name != canon_name:
                print("Header row {}: Name is {}, should be {}".format(idx, xl_name, canon_name))
        return False


def check_xl_bio_groups(xl_dataframe):
    bio_group_number = np.unique(xl_df.columns.get_level_values("Biological Group"))
    if bio_group_number == 2:
        return True
    else:
        print("Your Excel file contains {} biological groups, currently only 2 groups can be compared".format(bio_group_number))
        return False

                     
def check_normalisers(xl_dataframe, normaliser_filename):
    normalisers = get_candidate_normalisers(normaliser_filename)
    absent_normalisers = 0
    for normaliser in normalisers:
        if normaliser not in xl_dataframe_df.index:
            absent_normalisers += 1
            print("{} not found in Excel file".format(normaliser))
    if absent_normalisers > 0:
        return False
    else:
        return True

