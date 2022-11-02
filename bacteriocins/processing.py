## processing.py 

# All functions required to import, check and process a BIGSdb database export csv before data analysis and figure generation 

# Functions called within a Jupyter notebook 

# **Paper reference here** 

#####
##########__setup__##########
#####

## Importing libraries 

import sys
import os
import math
from pandas import Series, DataFrame
import pandas as pd
import numpy as np

from datetime import date

# setting the precision of calculations 
# to review - might not need this 
pd.set_option("precision", 3)

#####
##########__data_import__##########
#####

# data_import() finds the BIGSdb export and imports it 
# data must be in csv format 
# input_file variable must be defined in the notebook where process() is called 
# important that the read_csv function has the low_memory = False argument 
    # allows mixed data types in the allele designation columns (where there are "[S][I]" designations for example)
def data_import(input_file):
    print("Reading in file: " + input_file)
    return pd.read_csv(input_file, low_memory=False)

#####
##########__input_data_checks__##########
#####

# functions to check that the input data has the expected fields and that they are fully populated 

def csv_check(data):
    # defines some check functions, then calls them 
    # quits if a problem encountered 
    
    # spn_col_check checks that the expected columns are there 
    # specific for analyses with pneumococcal data only 
    def spn_col_check(data):
        # expected columns - required for analyses
        exp_columns = [
            "id",
            "final_serotype_designation",
            "clonal_complex",
            "year",
            "carriage_or_disease",
            "ST (MLST)",
            "source"
        ]
        # list of missing columns
        # added to if an expected column is not found 
        missing_columns = []
        for x in exp_columns:
            if x in data.columns:
                pass
            else:
                missing_columns.append(x)
        if len(missing_columns) == 0:
            # logging 
            print("Required columns found in df")
        else:
            print("Please check input csv, missing the following required columns: ")
            print(missing_columns)
            quit()
    # locus_check looks for curated locus columns 
    # defined as starting with "mj" 
    def locus_check(data):
        templocus = []
        for x in data.columns:
            if x.startswith("mj_"):
                templocus.append(x)
            else:
                pass
        # logging
        if len(templocus) == 0:
            print(
                "Please check input: no bacteriocin loci columns found (identified by mj_ prefix)"
                )
            print(
                "csv not processed"
            )
            quit()
        else:
            print("Found bacteriocin loci - continuing with processing")
    # this is not needed anymore 
    # leftover from when curation was incomplete 
    def uncurated_check(data):
        # all loci curated now - this list is empty 
        uncurated = []
        for x in uncurated:
            if x in data.columns:
                print(
                    "Found uncuratable bacteriocin loci in data export, do not export the following loci:"
                )
                print(uncurated)
                quit()
            else:
                pass

    spn_col_check(data)
    locus_check(data)
    #uncurated_check(data)
    
    # logging
    # if it gets to the end without quitting then the dataframe is fine and no action taken 
    print("Imported and checked data")


#####
##########__unique_clusters_and_loci__##########
#####

# useful functions callsed by downstream processing functions 

# function returning a list of all the loci columns in the dataframe 
def loci(data):

    # this function can be modified when the loci naming convention in BIGS is changed
    loci = [x for x in data.columns if x.startswith("mj")]
    return loci

# function returning a list of the unique bacteriocin clusters with loci in the given dataframe 

def unique_clusters(data):
    clusters = [x[3:6] for x in data.columns if x.startswith("mj")]
    uniques = []
    for x in clusters:
        if x in uniques:
            pass
        else:
            uniques.append(x)
    return uniques


#####
##########__processing_functions__##########
#####

#####__data_types__#####

# making sure that the allele designations are strings rather than integers
# need to be str format for manipulations in BacteriocinProfile class 

def df_types(data):
    for x in data.columns:
        if x.startswith("mj"):
            data["{}".format(x)] = data["{}".format(x)].astype(str)
        else:
            pass
    # logging 
    print("Adjusted data type of allele designations to string")


#####__column_names__#####

# adjusting col names to remove space in the ST column 

def col_names(data):
    
    data.rename(columns={"ST (MLST)": "mlst_st"}, inplace=True)
    
    # logging 
    print("Renamed columns")

#####__pre/post_vaccination__#####

# adds a column to the dataframe describing vaccine status 
# either in the pre or post PCV10 time period 
# defined based on the global vac_year variable 

def vac(data, vac_year):

    """Adds a column to the dataframe describing whether each isolate was collected in the pre or post vaccine time period
    Vac year is defined as a global variable in the setup cell at the start of the notebook"""

    vac_ser = Series(index=data.index, dtype="object")
    for x in data.index:
        if data["year"][x] != "" and data["year"][x] >= vac_year:
            vac_ser[x] = "post"
        elif data["year"][x] != "" and data["year"][x] < vac_year:
            vac_ser[x] = "pre"
        else:
            vac_ser[x] = "NA"
    # adds a column to the given dataframe describing vaccination time period 
    data["vaccination"] = vac_ser
    # logging 
    print("Defined pre/post vaccine periods based on year of sampling")

#####__disease_types__#####

# adds a column grouping diseases based on sample source 
# need a better breakdown than simple carriage/disease
# LRTI - from sputum or lung aspirate 
# OM - from middle ear fluid
# IPD - from sterile site (any source not listed above) 

def disease_breakdown(data):
    # function to return a list of the new disease code
    def disease_def(data):  
        templist = []
        for x in data.index:
            if data.carriage_or_disease[x] == "carriage":
                templist.append("Carriage")
            elif data.source[x] == "sputum":
                templist.append("LRT")
            elif data.source[x] == "lung aspirate":
                templist.append("LRT")
            elif data.source[x] == "middle ear fluid":
                templist.append("Otitis media")
            else:
                templist.append("Invasive")
        return templist
    # adds disease column to the given dataframe 
    data["disease_type"] = disease_def(data)
    # logging
    print("Defined disease type based on source field")

#####_locus_column_order__#####

# changing the order of the bacteriocin loci columns 
# to reflect the order the genes are found in in the bacteriocin clusters 
# otherwise it is alphabetical 
# this helps when describing partial or non-contiguous clusters 

# function returning the desired locus column order 
def adjust_order(data):
    # dict with the order of the non-alphabetic loci
    # this should be checked - done by hand
    # according to the order of loci in Ray's reference clusters 
    non_alph = {
        "sla": ["mj_slaA1", "mj_slaA2", "mj_slaA3", "mj_slaA4", "mj_slaA5",
            "mj_slaF", "mj_slaE", "mj_slaK", "mj_slaR", "mj_slaM", "mj_slaT",
        ],
        "slb": ["mj_slbF", "mj_slbG", "mj_slbE", "mj_slbA", "mj_slbM", "mj_slbT"],
        "slc": ["mj_slcA", "mj_slcX", "mj_slcL", "mj_slcT"],
        "sle": ["mj_sleM1", "mj_sleA1", "mj_sleA2", "mj_sleM2", "mj_sleM3",
            "mj_sleT", "mj_sleX1", "mj_sleF", "mj_sleG", "mj_sleX2"
        ],
        "slg": [ "mj_slgA1", "mj_slgA2", "mj_slgM", "mj_slgD",
            "mj_slgP1", "mj_slgT", "mj_slgP2",
        ],
        "slh": ["mj_slhP", "mj_slhR", "mj_slhK", "mj_slhF", "mj_slhE",
            "mj_slhG", "mj_slhX1", "mj_slhX2", "mj_slhA", "mj_slhB",
            "mj_slhT", "mj_slhC", "mj_slhI",
        ],
        "sli": [
            "mj_sliP","mj_sliR", "mj_sliK", "mj_sliF", "mj_sliE",
            "mj_sliG", "mj_sliA", "mj_sliB", "mj_sliT", "mj_sliC",
            "mj_sliI",
        ],
        "slj": ["mj_sljA1", "mj_sljL", "mj_sljP", "mj_sljT1", "mj_sljT2", 
                "mj_sljT3", "mj_sljA2"],
        "sls": ["mj_slsA", "mj_slsC", "mj_slsB1", "mj_slsB2", "mj_slsF",
            "mj_slsE", "mj_slsG", "mj_slsR", "mj_slsK",
        ],
        "ssa": ["mj_ssaA", "mj_ssaCD", "mj_ssaX1", "mj_ssaX2", "mj_ssaP", "mj_ssaX3"],
    }
    # specify the order of the loci
    locus_order = []
    for locus in loci(data):
        if locus[3:6] in non_alph.keys():
            if locus in locus_order:
                pass
            else:
                for entry in non_alph[locus[3:6]]:
                    locus_order.append(entry)
        else:
            locus_order.append(locus)

    # populating a list of all the column names of the df
    # first the non-locus names, then the locus names in the cluster order
    # as described in locus_order
    col_order = []
    for col_name in list(data.columns):
        if col_name.startswith("mj"):
            pass
        else:
            col_order.append(col_name)

    # add the locus order list onto the overall col order list
    col_order = col_order + locus_order
    return col_order

# a function for changing the order of the columns on the input dataframe
# calls adjust_order(data) 
def col_order(data):
    data = data[adjust_order(data)]
    print(
        "Adjusted locus order in curated data export to reflect order found in reference clusters"
    )
    return data


#####__BacteriocinProfile_class__#####

# a class to hold information about the loci of a cluster
# could be more sophisticated - is better in Manatee 
# not worrying about it now, it works  
class bacteriocin_profile:
    def __init__(self, dict):

        self.cluster = list(dict.keys())[0][:3]

        self.allelic_profile = "-".join(dict.values())

        self.partfull_profile = "-".join(
            ["/" if v == "0" else k[3:] for k, v in dict.items()]
        )

        self.category = cluster_category(self.partfull_profile)

        # hashed out the status attribute 
        # became more complex after the incorporation of cluster contiguity checks 
        # self.status = True if self.category in ["P", "Fl"] else False

# allele_dict generates a dictionary with loci and alleles for a given cluster from the dataframe of curated data
# for input into the bacteriocin_profile class
def allele_dict(data, cluster, index):

    loci = [x[3:] for x in data.columns if cluster in x]
    output = {}
    for col in loci:
        output[col] = data["mj_{}".format(col)][index]

    return output

# function for generating the generalised full/partial status attribute in the bacteriocin_profile class
def cluster_category(profile):
    # profile arg is a partfull profile e.g. /-B-C for a streptococcin

    if len(set(profile.split("-"))) == 1:
        return "A"
    elif "/" in set(profile.split("-")):
        if len(set(profile.split("-"))) == 2:
            return "Fg"
        else:
            return "P"
    else:
        return "Fl"

# could write some exceptions in here if needed
# this follows the simple rule that if a single locus is observed then that counts as a fragment, not as a partial
 

#####__BacteriocinProfile_columns__#####

# functions using the BacteriocinProfile class to generate new columns in the dataframe 

# this returns a dataframe with the index of the dataframe with columns for all the bacteriocin_profile attributes
# this can be merged onto the dataframe
def cluster_cols(data, cluster):

    cluster_df = DataFrame(
        index=data.index,
        columns=[
            "{}_a_profile".format(cluster),
            "{}_profile".format(cluster),
            "{}_category".format(cluster),
            #"{}_status".format(cluster),
        ],
    )

    for x in data.index:
        cluster_df["{}_a_profile".format(cluster)][x] = bacteriocin_profile(
            allele_dict(data, cluster, x)
        ).allelic_profile
        cluster_df["{}_profile".format(cluster)][x] = bacteriocin_profile(
            allele_dict(data, cluster, x)
        ).partfull_profile
        cluster_df["{}_category".format(cluster)][x] = bacteriocin_profile(
            allele_dict(data, cluster, x)
        ).category
        # cluster_df["{}_status".format(cluster)][x] = bacteriocin_profile(
        #   allele_dict(data, cluster, x)
        # ).status
    
    # logging 
    print("Generated profile and category columns for {}".format(cluster))
    
    return cluster_df

# this function calls cluster_cols for all the bacteriocins of interest and merges the new cols onto the given dataframe
# returns a separate object
def profile_cols(data):

    output = data
    for x in unique_clusters(data):
        output = output.merge(
            cluster_cols(data, x), how="outer", left_index=True, right_index=True
        )
    
    # logging 
    print("Added profile columns for bacteriocin clusters to the dataframe")
    
    return output


#####__Bacteriocin_status_columns__#####

# Adding a column for each bacteriocin describing the status of the bacteriocin cluster 
# Boolean type - True if the bacteriocin is present and should be counted in analysis
# False if no loci of the bacteriocin were detected, if the cluster is a fragment, or if the cluster is non-contiguous 

# Uses the contiguity cat script output to find clusters which are not contiguous and should not be counted 
# Importing the contiguity cat output from the local directory 
def contiguity_cat_import():
    
    # importing the longform output from contiguity_cat.py 
    # with 2500bp thresholds for adjacent loci and for EOCs
    I_cat = pd.read_csv("../../data/contiguity_cat_outputs/cluster_cont_VICE_2500_2500.csv")
    K_cat = pd.read_csv("../../data/contiguity_cat_outputs/cluster_cont_Kenya_2500_2500.csv")

    # joining them together 
    cluster_cat = pd.concat([I_cat, K_cat])
    cluster_cat.set_index("id", inplace = True)
    return cluster_cat

# cluster status returns a Series object with the status information 
def cluster_status(data, cluster, cat_df):
    
    status_col = Series(index = data.index, dtype = bool)
    for row in data.index:
        isolate = data.id[row]
        if cat_df["{}_contiguitycat".format(cluster)][isolate].startswith("Non-contiguous"):
            status_col[row] = False
        elif cat_df["{}_contiguitycat".format(cluster)][isolate] == "Absent":
            status_col[row] = False
        else:
            status_col[row] = True
    
    return status_col 

# adds a column for the T/F status of each cluster to the dataframe
# adds it directly to the given dataframe, does not return a new dataframe 
def status_cols(data, cat_df):
    
    for cluster in unique_clusters(data):
        data["{}_status".format(cluster)] = cluster_status(data, cluster, cat_df)
    print("Added T/F status columns for all bacteriocins in the dataset")


#####___Bacteriocin_cluster_count___#####

# Counting how many different bacteriocin clusters are detectable in each genome 
# based on the status columns - how many clusters are True?
# must be called after the dataframe has been given its status columns 
# one new column to the given dataframe with the count 

def cluster_count(data):

    # generating a list of the status columns
    status_cols = []
    for x in data.columns:
        if x.endswith("status"):
            status_cols.append(x)

    # slice of the dataframe with only status columns
    statusdf = data[status_cols]
    count_col = Series(index=data.index, dtype="object")

    for x in statusdf.index:
        temp = statusdf.iloc[x].value_counts()
        if True in temp.index:
            count_col[x] = temp[True]
        else:
            count_col[x] = 0
    count_col_int = count_col.astype(int)
    data["cluster_count"] = count_col_int
    # logging 
    print("Generated cluster count column")

#####___Bacteriocin_repertoire___#####

# Related to the cluster count function - what combination of clusters is present in each genome?
# Generates a long string describing the combination of bacteriocins 
# Not intended to be human readable, but can be plotted as heatmaps 
# adds a new column to the given dataframe 

def bact_rep(data):

    # this function generates the repertoire for a single genome 
    def repertoire(isolate, data):

        clusters = []
        for x in unique_clusters(data):
            if data["{}_status".format(x)][isolate] == True:
                clusters.append(x)
            else:
                clusters.append("/")

        rep = "-".join(clusters)
        return rep

    # calling repertoire() for all rows in the dataframe
    repertoire_col = Series(index=data.index, dtype="object")

    for row in data.index:
        repertoire_col[row] = repertoire(row, data)
    data["bacteriocin_repertoire"] = repertoire_col
    # logging 
    print("Generated bacteriocin repertoire column")


#####
##########__process()__##########
#####

# A single function that imports and processes a dataframe of genomic data with bacteriocin gene annotations 
# returns the dataframe with all the additional columns added ready for analysis section 

def process(input_file, vac_year):
    # imports the dataframe
    rawdf = data_import(input_file)

    # checking the columns and adjusting types/names
    csv_check(rawdf)
    df_types(rawdf)
    col_names(rawdf)

    # adding columns with processed vaccine/disease status
    vac(rawdf, vac_year)
    disease_breakdown(rawdf)
    
    # adjusting the order of the loci columns to reflect locus order in bacteriocin clusters
    rawdf = col_order(rawdf)
    
    # adding columns with bacteriocin profiles using the profile class attributes
    rawdf = profile_cols(rawdf)
    
    # adding a status column for each bacteriocin 
    # incorporates contiguity cat outputs 
    status_cols(rawdf, contiguity_cat_import())
    
    # cluster count - how many bacteriocins are detectable in each genome?
    cluster_count(rawdf)

    # bact_rep(df) adds a column describing the combination of bacteriocins which are detectable in the genome
    bact_rep(rawdf)

    print("Processed data")
    return rawdf


