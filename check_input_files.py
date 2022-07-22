import pandas as pd
import numpy as np
import sys

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
        print("Your Excel file contains {} biological groups, currently only 2 groups can be compared.".format(bio_group_number))
        return False

    
def normaliser_warning(list_of_normalisers):
    if len(list_of_normalisers) > 11:
        continue_ask = input("This comparison will take more than 40 minutes, and your system may run out of memory.\nWould you like to continue? (y/n) ")
        if continue_ask == "y":
            pass
        else:
            sys.exit()

            
def get_candidate_normalisers(filename, sheet_name="normalisers"):
    normaliser_df = pd.read_excel(filename, sheet_name=sheet_name)
    return normaliser_df


def get_normaliser_locations(normaliser_df):
    normaliser_locations = []
    for normaliser in normalisers:
        if normaliser not in pcr_dataframe.index:
            
            print("{} not found in Excel file".format(normaliser))
        else:
            normaliser_locations.append(np.where(pcr_dataframe.index == normaliser)[0])
    return normaliser_locations

    
def check_normalisers(pcr_dataframe, normaliser_filename):
    normalisers = get_candidate_normalisers(normaliser_filename)
    normaliser_warning(normalisers)
    normaliser_index_locations = get_normaliser_locations(normalisers)
    if len(normaliser_index_locations) == len(normalisers):
        return True
    else:
        return False


def all_checks(xl_dataframe, normaliser_filename):
    checks = (check_xl_column_names(xl_dataframe), 
              check_xl_bio_groups(xl_dataframe),
              check_normalisers(xl_dataframe, normaliser_filename))
    if np.all(checks):
        return True
    else:
        return False
    
