## analysis.py 

# Functions for further processing, summarisation and visualisation of a processed BIGSdb dataset export

# File organised according to figures, tables and supplementary materials generated for the paper (reference here when it exists)
# Assumes that the dataframe fed to each function has been imported and processed using process.py first 

###
#####_____SETUP_____#####
###

## Importing libraries used by functions defined here 
import sys
import os
import math
from this import d
from pandas import Series, DataFrame
import pandas as pd
import numpy as np

from datetime import date

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import seaborn as sns

## Setting up colour palettes 

# use the seaborn colorblind palette as default 
sns.set_palette("colorblind")

# defining a palette for general use in plotting functions where bars are broken down by disease type 
disease_pal = ["#38AECC", "#db5461", "#c0da74", "#F4A259"]

#####
#####           
##########____________________GENERAL_FUNCTIONS____________________##########
#####
#####

## unique_clusters() and unique_loci()
# generates a list of the clusters and loci included in the dataset export 
# uses the prefic "mj_" to recognise bacteriocin loci 

def unique_clusters(data):
    clusters = [x[3:6] for x in data.columns if x.startswith("mj")]
    uniques = []
    for x in clusters:
        if x in uniques:
            pass
        else:
            uniques.append(x)
    return uniques

def loci(data):

    # this function can be modified when the loci naming convention in BIGS is changed
    loci = [x for x in data.columns if x.startswith("mj")]
    return loci

## present_clusters() and absent_clusters()
# generating lists of clusters which were detectable in a given dataset 

def present_clusters(data):
    output = [x for x in unique_clusters(data) if data["{}_status".format(x)].any()]
    return output


def absent_clusters(data):
    output = [x for x in unique_clusters(data) if not data["{}_status".format(x)].any()]
    return output

## general chi-square function 
# for the chi-square test to compare two proportions 
# used by various functions to detect differences in bacteriocim prevalence between subsets of the dataset 

def chisq(proc_data):

    # data is a dataframe, returned by a processing function, indexed by bacteriocin
    # the processing function is different for every application of this function - defined as required below 
    # frequencies of the 20 bacteriocins in "condition 1" and "condition 2", could be e.g. pre/post vac
    # return a dataframe with each bacteriocin, saying either "significant", "not significant" or "sample size"
    # distinguishes bacteriocins which are not significantly different from cases where the sample size of the bacteriocin is not sufficient to assess significance 

    output = DataFrame(index=proc_data.index, columns=["chisq_result", "chisq_value"])

    for x in proc_data.index:
        O1 = proc_data.O1[x]
        O2 = proc_data.O2[x]
        O3 = proc_data.O3[x]
        O4 = proc_data.O4[x]

        RM1 = O1 + O2
        RM2 = O3 + O4
        CM1 = O1 + O3
        CM2 = O2 + O4
        N = RM1 + RM2

        E1 = (RM1 * CM1) / N
        E2 = (RM1 * CM2) / N
        E3 = (RM2 * CM1) / N
        E4 = (RM2 * CM2) / N

        # if any value of E is <5, the sample size is too low 
        if E1 < 5:
            output.chisq_result[x] = "LowN"
            continue
        elif E2 < 5:
            output.chisq_result[x] = "LowN"
            continue
        elif E3 < 5:
            output.chisq_result[x] = "LowN"
            continue
        elif E4 < 5:
            output.chisq_result[x] = "LowN"
            continue
        else:
            pass

        # Yates correction continuity factor 
        # (O1 - E1 - 0.5)**2/E1
        X1 = ((O1 - E1) ** 2) / E1
        X2 = ((O2 - E2) ** 2) / E2
        X3 = ((O3 - E3) ** 2) / E3
        X4 = ((O4 - E4) ** 2) / E4

        overallX = X1 + X2 + X3 + X4

        if overallX > 10.83:
            output.chisq_result[x] = "P < 0.001"
        elif overallX > 6.64:
            output.chisq_result[x] = "P < 0.01"
        elif overallX > 3.84:
            output.chisq_result[x] = "P < 0.05"
        else:
            output.chisq_result[x] = "NS"

        output.chisq_value[x] = overallX
        print("Chisq N = " + str(N))

    return output


#####
#####           
##########____________________DATASET_SUMMARIES____________________##########
#####
#####

# Generating tables and figures that summarise the datasets 
# Not handling the bacteriocin loci yet 

###
#####___TABLE_1_____#####
###

# summarising lots of fields in the processed datasets
# including metadata and lineage fields 

def summary_table(data):

    subdata = [data[data.country == "Iceland"], data[data.country == "Kenya"]]

    summaries = []
    for subdf in subdata:

        total = len(subdf)

        # breakdown by carriage/disease
        carriage_count = len(subdf[subdf.carriage_or_disease == "carriage"])
        disease_count = len(subdf[subdf.carriage_or_disease == "disease"])
        invasive_count = len(subdf[subdf.disease_type == "Invasive"])

        lrti_count = 0
        if "LRT" in subdf.disease_type.unique():
            lrti_count = len(subdf[subdf.disease_type == "LRT"])

        om_count = 0
        if "Otitis media" in subdf.disease_type.unique():
            om_count = len(subdf[subdf.disease_type == "Otitis media"])

        # breakdown by pre/post vaccination period
        prevac_count = len(subdf[subdf.vaccination == "pre"])
        postvac_count = len(subdf[subdf.vaccination == "post"])

        # breakdown by gender of patient
        male_count = len(subdf[subdf.gender == "male"])
        female_count = len(subdf[subdf.gender == "female"])

        # number of unique clonal complexes
        unique_cc = len(subdf.clonal_complex.unique())

        # unique clonal complexes excluding singletons
        sing = [x for x in subdf.clonal_complex.unique() if x.startswith("Sing")]
        unique_cc_exclsing = unique_cc - len(sing)

        # unique STs
        unique_st = len(subdf.mlst_st.value_counts())

        # unique serotypes
        unique_serotype = len(subdf.final_serotype_designation.value_counts())

        values = [
            total,
            carriage_count,
            disease_count,
            invasive_count,
            lrti_count,
            om_count,
            prevac_count,
            postvac_count,
            male_count,
            female_count,
            unique_cc,
            unique_cc_exclsing,
            unique_st,
            unique_serotype,
        ]
        index = [
            "Total",
            "Carriage",
            "Disease",
            "Invasive",
            "LRTI",
            "Otitis media",
            "Pre-vaccination",
            "Post-vaccination",
            "Male",
            "Female",
            "Unique clonal complexes",
            "(excluding singletons)",
            "Unique STs",
            "Unique serotypes",
        ]

        summary = Series(values, index=index)
        summary = DataFrame(summary)
        if list(subdf.country.unique()) == ["Iceland"]:
            summary.rename(columns={0: "Iceland"}, inplace=True)
        elif list(subdf.country.unique()) == ["Kenya"]:
            summary.rename(columns={0: "Kenya"}, inplace=True)
        summaries.append(summary)

    output = pd.concat(summaries, axis=1)

    return output

###
#####___SHARED_CCs_AND_STs_____#####
###

# functions which return the number of CCs or STs which are found in both datasets 

def shared_CCs(data):
    shared_list = []

    I_unique = data[data.country == "Iceland"].clonal_complex.unique()
    K_unique = data[data.country == "Kenya"].clonal_complex.unique()
    print("Iceland:")
    print(len(I_unique))
    print("Kenya:")
    print(len(K_unique))

    for x in I_unique:
        if x in K_unique:
            shared_list.append(x)
        else:
            pass

    # returns the length of the list i.e. how many ccs are found in both datasets
    return len(shared_list)

    # this would return the list of which ccs are found in both datasets
    # return(shared_list)

def shared_STs(data):
    shared_list = []

    I_unique = data[data.country == "Iceland"].mlst_st.unique()
    K_unique = data[data.country == "Kenya"].mlst_st.unique()
    print("Iceland:")
    print(len(I_unique))
    print("Kenya:")
    print(len(K_unique))

    for x in I_unique:
        if x in K_unique:
            shared_list.append(x)
        else:
            pass

    return len(shared_list)


###
#####___TABLE_2_____#####
###

# generates a table of the most common CCs in the given dataset 
# top X most common determined by threshold argument 
# any not in top X pooled to "other"
# can be run on the whole dataset, or on subsets
# most useful when run on one country at a time 

def cc_freq_proc(data, threshold):
    # value counts of the clonal complexes in the dataset 
    vc = data.clonal_complex.value_counts()
    # subset according to the threshold value 
    top_vc = vc.head(threshold)
    
    # get list of ccs 
    ccs = list(top_vc.index)
    ccs.append("Other CCs")
    ccs.append("Singletons")
    
    
    output = DataFrame(columns = ["CC", "n (%)"])
    output["CC"] = ccs
    
    total = len(data)
    
    # gerenating the other/singleton values
    singleton_count = 0
    other_count = 0
    for row in vc.index:
        if row in top_vc.index:
            pass
        elif row.startswith("Sing"):
            singleton_count = singleton_count + vc[row]
        else:
            other_count = other_count + vc[row]
    singleton_prop = (singleton_count/total)*100
    other_prop = (other_count/total)*100
    
    for row in output.index:
        
        cc = output["CC"][row]
        
        if cc == "Other CCs":
            output["n (%)"][row] = str(other_count) + " (" + str(round(other_prop, 1)) + "%)"
            #output.Frequency[row] = other_count
            #output["Proportion of dataset"][row] = other_prop
        elif cc == "Singletons":
            output["n (%)"][row] = str(singleton_count) + " (" + str(round(singleton_prop, 1)) + "%)"
        else:
            count = vc[cc] 
            prop = (count/total)*100
            output["n (%)"][row] = str(count) + " (" + str(round(prop, 1)) + "%)"
            
            #output.Frequency[row] = count
            #output["Proportion of dataset"][row] = prop 
        
    return output


###
#####___TABLE_3_____#####
###

def serotype_freq_proc(data, threshold):
    # value counts of the clonal complexes in the dataset 
    vc = data.final_serotype_designation.value_counts()
    # subset according to the threshold value 
    top_vc = vc.head(threshold)
    
    # get list of ccs 
    serotypes = list(top_vc.index)
    serotypes.append("Other serotypes")
    
    output = DataFrame(columns = ["Serotype", "n (%)"])
    output["Serotype"] = serotypes
    
    total = len(data)
    
    # gerenating the other values
    other_count = 0
    for row in vc.index:
        if row in top_vc.index:
            pass
        else:
            other_count = other_count + vc[row]
    other_prop = (other_count/total)*100
    
    for row in output.index:
        
        serotype = output["Serotype"][row]
        
        if serotype == "Other serotypes":
            output["n (%)"][row] = str(other_count) + " (" + str(round(other_prop, 1)) + "%)"
        else:
            count = vc[serotype] 
            prop = (count/total)*100
            output["n (%)"][row] = str(count) + " (" + str(round(prop, 1)) + "%)"
        
    return output

###
#####___FIGURE_1A_____#####
###

# Figure showing the distribution of sampling year 
# 1 panel per dataset, both returned as a single figure 

## disease_legend
# function to generate a legend for any stacked bar plots showing the colours representing each disease type 
# used in figure 1A and 1B

def dis_cat_leg(ax, coords):
    
    # different rectangles for the various disease categories 
    carbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[0], edgecolor="none")
    invbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[1], edgecolor="none")
    lrtbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[2], edgecolor="none")
    ombar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[3], edgecolor="none")
    
    # draws the legend at the specified coordinates 
     # on the specified axis 
    ax.legend(
        [carbar, invbar, lrtbar, ombar],
        ["Carriage", "IPD", "LRTI", "OM"],
        loc="upper right",
        ncol=4,
        prop={"size": 10},
        bbox_to_anchor=coords,
    )

# first a processing function:
    # Reasonably complex as generating a stacked bar chart 

def group_disease(data):
    # this is a non-melted type dataframe
     # one column for each disease category
     # consecutively adding each category to the last to get the stacking

    car_stack = data[data.disease_type == "Carriage"].year.value_counts()

    inv_stack = car_stack.add(
        data[data.disease_type == "Invasive"].year.value_counts(), fill_value=0
    )

    grouped_df = DataFrame({"car_stack": car_stack, "inv_stack": inv_stack})

    if "LRT" in data.disease_type.unique():
        lrt_stack = inv_stack.add(
            data[data.disease_type == "LRT"].year.value_counts(), fill_value=0
        )
        grouped_df["lrt_stack"] = lrt_stack

        om_stack = lrt_stack.add(
            data[data.disease_type == "Otitis media"].year.value_counts(), fill_value=0
        )
        grouped_df["om_stack"] = om_stack
    else:
        pass
    grouped_df.reset_index(inplace=True)
    grouped_df.rename(columns={"index": "year"}, inplace=True)
    return grouped_df

# then the plotting function:

def stack_year_plot(data):

    # splitting out I_df and K_df from the overall dataset 
    I_df = data[data.country == "Iceland"]
    K_df = data[data.country == "Kenya"]

    # Aesthetics
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    # Setting up axes and subplots
    stackfig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 5))
    plt.subplots_adjust(hspace=0.4)

    # Plotting different graphs onto the axes to achieve stacked bars
    # got here using group_disease(data)
    sns.barplot(
        x="year", y="om_stack", data=group_disease(I_df), color=disease_pal[3], ax=ax1
    )
    sns.barplot(
        x="year", y="lrt_stack", data=group_disease(I_df), color=disease_pal[2], ax=ax1
    )
    sns.barplot(
        x="year", y="inv_stack", data=group_disease(I_df), color=disease_pal[1], ax=ax1
    )
    sns.barplot(
        x="year", y="car_stack", data=group_disease(I_df), color=disease_pal[0], ax=ax1
    )

    sns.barplot(
        x="year", y="inv_stack", data=group_disease(K_df), color=disease_pal[1], ax=ax2
    )
    sns.barplot(
        x="year", y="car_stack", data=group_disease(K_df), color=disease_pal[0], ax=ax2
    )

    # Manually defining a legend using dis_cat_leg()
    leg1 = dis_cat_leg(ax1, (0.88, -1.7))

    # Titles
    stackfig.axes[0].set_title("Iceland", fontsize=12)
    stackfig.axes[1].set_title("Kenya", fontsize=12)

    # X axis labels
    stackfig.axes[0].set_xlabel("", fontsize=12)
    stackfig.axes[1].set_xlabel("Year", fontsize=12)

    # Y axis labels
    for x in stackfig.axes:
        x.set_ylabel("Pneumococci (n)", fontsize=12)

    # adding lines for vac year
    ax1.axvline(2.5, 0, 1, linestyle="--", color="black")
    ax2.axvline(8.5, 0, 1, linestyle="--", color="black")

    return stackfig


###
#####___FIGURE_1B_____#####
###

# a figure showing the distribution of patient age in the two datasets 
# again stacked by disease type and using the same figure legend as above

# processing function:
# generating data for a stacked bar, so complicated again 
 # as for other stacked bars, each layer is generated by adding consecutive y values from each category 
# also bins the ages 
def age_stack_process(data):

    # bins with 0-5, 5-10 and then 10 year bins up to 100
    bins = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # setting up the bins on the output dataframe
    outindex = (pd.cut(data.age_yrs, bins)).value_counts(sort=False).index
    output = DataFrame(index=outindex)
    output.reset_index(inplace=True)
    output.columns = ["age"]

    # the different stack categories
    stacks = data.disease_type.unique()
    for x in stacks:
        tempdf = data[data.disease_type == x]
        tempcut = pd.cut(tempdf.age_yrs, bins)
        tempcount = DataFrame(tempcut.value_counts(sort=False)).reset_index()
        tempcount.columns = ["age", "{}_count".format(x)]

        output = pd.merge(output, tempcount, on="age")

    # carriage is the first of the stack categories, not added to anything
    output.rename(columns={"Carriage_count": "Carriage"}, inplace=True)

    # generating the invasive category by adding it to carriage
    inv_stack = []
    for x in output.index:
        inv_stack.append((output.Carriage[x]) + (output.Invasive_count[x]))
    output["Invasive Disease"] = inv_stack

    # alternate data prep for the two datasets
    # reflects the different disease sampling in Iceland and Kenya
    if list(data.country.unique()) == ["Iceland"]:
        lrt_stack = []
        om_stack = []
        for x in output.index:
            lrt_stack.append((output["Invasive Disease"][x]) + (output.LRT_count[x]))
            om_stack.append(
                (output["Invasive Disease"][x])
                + (output.LRT_count[x])
                + output["Otitis media_count"][x]
            )
        output["LRT Infection"] = lrt_stack
        output["Otitis Media"] = om_stack

    elif list(data.country.unique()) == ["Kenya"]:
        output["LRT Infection"] = len(output) * [0]
        output["Otitis Media"] = len(output) * [0]

    # dropping columns from output that are not needed
    bin_columns = [x for x in output.columns if x.endswith("count")]
    output.drop(columns=bin_columns, inplace=True)

    # generating the final output in a convenient format for plotting (melted)
    varlist = ["Carriage", "Invasive Disease", "LRT Infection", "Otitis Media"]
    output_m = pd.melt(
        output, id_vars=["age"], value_vars=[x for x in output.columns if x in varlist]
    )

    output_m["Dataset"] = (len(output_m)) * list(data.country.unique())
    return output_m

# then the plotting functions 

# generalised figure for plotting a single subplot
 # plots onto the specified axis using the argument data (always the output of age_stack_process(data))
def age_subplot(ax, data): 
    sns.barplot(
        x="age",
        y="value",
        data=data,
        hue="variable",
        palette=[disease_pal[3], disease_pal[2], disease_pal[1], disease_pal[0]],
        hue_order=["Otitis Media", "LRT Infection", "Invasive Disease", "Carriage"],
        dodge=False,
        ax=ax,
    )
# function plotting the subplots onto a set of axes 
    # complex because of smaller zoomed in subplots on the main panels showing older (less well-sampled) age groups 
def age_plot(data):

    # set up
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    # setting up axes
    agefig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 5))
    plt.subplots_adjust(hspace=0.55)
    axes = [ax1, ax2]

    # pulling out I_df and K_df
    I_df = data[data.country == "Iceland"]
    K_df = data[data.country == "Kenya"]

    # Plotting different graphs onto the axes to achieve stacked bars
    # Icelandic plot
    age_subplot(ax1, age_stack_process(I_df))
    
    # Kenyan plot
    age_subplot(ax2, age_stack_process(K_df))
    
    # aesthetics
    axes[0].set_title("Iceland", fontsize=12)
    axes[1].set_title("Kenya", fontsize=12)

    for x in axes:
        x.set(xlim=[10.5, -0.5])
        x.set_ylabel("Pneumococci (n)", fontsize=12)
        x.spines["right"].set_visible(True)
        x.spines["top"].set_visible(True)
        x.invert_xaxis()
        x.set_xticklabels(
            [
                "< 5",
                "5-9.9",
                "10-19.9",
                "20-29.9",
                "30-39.9",
                "40-49.9",
                "50-59.9",
                "60-69.9",
                "70-79.9",
                "80-89.9",
                "> 90",
            ],
            fontsize=10,
            rotation=45,
        )
        x.legend_.remove()
        x.yaxis.set_major_locator(plt.MaxNLocator(5))

    ax1.set_xlabel("", fontsize=1)
    ax2.set_xlabel("Age (years)", fontsize=12)

    ax1.set(ylim=[0, 1300])
    ax2.set(ylim=[0, 2100])

    # Manually defining a legend using def_cat_leg() 
    leg1 = dis_cat_leg(axes[1], (1, -0.6))

    # Adding zoomed in sub-plots
    I_subset = (age_stack_process(I_df).age.unique())[2:]
    K_subset = (age_stack_process(K_df).age.unique())[3:]

    ax3 = plt.axes([0.349, 0.75, 0.55, 0.13])  # this is [xcoord, ycoord, xsize, ysize]
    ax4 = plt.axes([0.41, 0.291, 0.49, 0.13])

    # plotting Icelandic subplot
    age_subplot(ax3, age_stack_process(I_df)[age_stack_process(I_df).age.isin(I_subset)])

    # plotting Kenyan subplot    
    age_subplot(ax4, age_stack_process(K_df)[age_stack_process(K_df).age.isin(K_subset)])
    
    # subplot aesthetics 
    subplots = [ax3, ax4]
    for ax in subplots:
        ax.invert_xaxis()
        ax.legend_.remove()
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))
     
    # subplots - x axis labels 
    ax3.set(xlim=[1.5, 10.5], ylim=[0, 130], ylabel="", xlabel="")
    ax3.set_xticklabels(
        [
            "",
            "",
            "10-19.9",
            "20-29.9",
            "30-39.9",
            "40-49.9",
            "50-59.9",
            "60-69.9",
            "70-79.9",
            "80-89.9",
            "> 90",
        ],
        fontsize=9,
        rotation=45,
    )
    
    ax4.set(xlim=[2.5, 10.5], ylim=[0, 130], ylabel="", xlabel="")
    ax4.set_xticklabels(
        [
            "",
            "",
            "",
            "20-29.9",
            "30-39.9",
            "40-49.9",
            "50-59.9",
            "60-69.9",
            "70-79.9",
            "80-89.9",
            "> 90",
        ],
        fontsize=9,
        rotation=45,
    )    

    return agefig


###
#####___Serotype_distribution_____#####
###

# thesis only: code for a stacked bar chart showing the distribution of serotypes 




#####
#####           
##########____________________BACTERIOCIN_PREVALENCE____________________##########
#####
#####

# looking at the overall frequency of full and partial bacteriocin clusters in the two datasets 
# including recognising any bacteriocins that are significantly more or less common in one dataset using the chi-square test 
# general chi-square function defined in the general functions section 


###
#####_____PROCESSING/GENERAL_FUNCTIONS_____#####
###

## cluster summary function
# dataset processing to calculate the prevalence of each cluster in the dataset
# also to organise the data in a useful format for plotting a grouped, stacked bar chart 
def cluster_summary(data):
    # data argument is either Icelandic or Kenyan subset of main df 
    clustersum_dictlist = []

    def cluster_prev(ref, data):
        # ref is a three letter cluster reference from the unique_clusters() general function
        nonlocal clustersum_dictlist
        tempstatus = data["{}_status".format(ref)]
        temp_partial = data["{}_category".format(ref)]

        a = "oops"

        if True in tempstatus.value_counts().index:
            a = tempstatus.value_counts()[True]
        else:
            a = 0

        total = len(data)
        b = round(((a / total) * 100), 1)

        c = "oops"

        if "P" in temp_partial.value_counts().index:
            c = temp_partial.value_counts()["P"]
        else:
            c = 0

        d = round(((c / total) * 100), 1)

        tempdict = {
            "cluster": ref,
            "prev_count": a,
            "prev_perc": b,
            "part_count": c,
            "part_perc": d,
        }
        clustersum_dictlist.append(tempdict)

    for x in unique_clusters(
        data
    ):  # run cluster_prev for every cluster in the dataset
        cluster_prev(x, data)

    return DataFrame(clustersum_dictlist)

## chi-square processing function 
# uses cluster_summary() to generate prevalences 
def chisq_proc_dataset(data):
    # function preparing data for input into general chisq function

    I_df = data[data.country == "Iceland"]
    K_df = data[data.country == "Kenya"]

    output = DataFrame(index=unique_clusters(data))

    output["O1"] = cluster_summary(I_df).set_index("cluster").prev_count
    output["O2"] = cluster_summary(K_df).set_index("cluster").prev_count

    O3 = Series(index=output.index, dtype=int)
    O4 = Series(index=output.index, dtype=int)
    for x in output.index:
        O3[x] = len(I_df) - output.O1[x]
        O4[x] = len(K_df) - output.O2[x]

    output["O3"] = O3
    output["O4"] = O4

    return output


###
#####___SUPPLEMENTARY_TABLE_5_____#####
###

# function generating the table 
# classes clusters as "full", "partial" or "fragment" according to their category from BacteriocinProfile (dependant on how many of the expected loci were detected)
def supp_table_5_6(data):
    def single_pivot(cluster, data):

        temp_pivot = pd.pivot_table(
            data,
            values="id",
            index=[
                "country",
                "{}_profile".format(cluster),
                "{}_category".format(cluster),
            ],
            aggfunc=len,
        )
        temp_pivot["Bacteriocin"] = len(temp_pivot) * ["{}".format(cluster)]
        temp_pivot.reset_index(inplace=True)
        # rename columns here
        temp_pivot.rename(
            columns={
                "country": "Country",
                "{}_profile".format(cluster): "Profile",
                "{}_category".format(cluster): "Category",
                "id": "Frequency",
            },
            inplace=True,
        )
        temp_pivot = temp_pivot[
            ["Bacteriocin", "Country", "Profile", "Category",  "Frequency"]
        ]
        temp_pivot = temp_pivot.reindex(index=temp_pivot.index[::-1])
        
        # removing all the absent clusters
        # retaining the non-contiguous clusters, they are dropped at a later point 
        indexNames = temp_pivot[temp_pivot['Category'] == "A"].index
        temp_pivot.drop(indexNames, inplace=True)
        
        temp_pivot.replace("Fl", "Full", inplace = True)
        temp_pivot.replace("P", "Partial", inplace = True)
        temp_pivot.replace("Fg", "Fragment", inplace = True)
        
        return temp_pivot

    pivot_list = []
    for cluster in present_clusters(data):
        pivot_list.append(single_pivot(cluster, data))

    return pd.concat(pivot_list, ignore_index=True)

###
#####_____FIGURE_2A_____#####
###

## functions to generate a grouped stacked bar chart of the overall prevalence of each bacteriocin in each dataset
# including finding bacteriocins which are significantly more/less common in one dataset 

## plotting function 
def grouped_prev_bar(data):
    # style settings
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    data1 = cluster_summary(data[data.country == "Iceland"])
    data2 = cluster_summary(data[data.country == "Kenya"])

    # function to combine the data from the two locations into one dataframe for plotting 
    def tidy_prev_data(data1, data2):
        data1["Dataset"] = 20 * ["Iceland"]
        data2["Dataset"] = 20 * ["Kenya"]

        merged_data = pd.concat([data1, data2])
        merged_data = merged_data.drop(["prev_count"], axis=1)

        return merged_data
    # tidying up the data for use in plotting
    plot_data = tidy_prev_data(data1, data2)

    # plotting the bar chart using seaborn catplot to acheive grouped bars
    prev_plot = sns.catplot(
        x="cluster",
        y="prev_perc",
        hue="Dataset",
        data=plot_data,
        kind="bar",
        height=4,
        aspect=2.5,
        legend_out=False,
    )

    # overlaying bars showing the proportion of isolates with partial clusters of each bacteriocin
    sns.barplot(
        x="cluster",
        y="part_perc",
        hue="Dataset",
        palette=["#9dc7c8", "#fbbb4b"],
        data=plot_data,
        ax=prev_plot.ax,
    )

    # removing the legend
    # legend added manually in Affinity Designer 
    prev_plot.ax.legend_.remove()
    
    # labelling axes and ticks
    prev_plot.ax.set_ylabel("Prevalence (%)", fontsize=12)
    prev_plot.ax.set_xlabel("")
    bact_names = [
        "Cib",
        "Streptococcin A",
        "Streptococcin B",
        "Streptococcin C",
        "Streptococcin D",
        "Streptococcin E",
        "Streptocyclicin",
        "Streptolancidin A",
        "Streptolancidin B",
        "Streptolancidin C",
        "Streptolancidin D",
        "Streptolancidin E",
        "Streptolancidin F",
        "Streptolancidin G",
        "Streptolancidin H",
        "Streptolancidin I",
        "Streptolancidin J",
        "Streptolancidin K",
        "Streptolassin",
        "Streptosactin",
    ]
    prev_plot.ax.set_xticklabels(
        bact_names, fontsize=12, rotation=45, ha="right", rotation_mode="anchor"
    )
    prev_plot.ax.set_yticklabels(prev_plot.ax.get_yticklabels(), fontsize=12)

    # box around plot
    prev_plot.ax.spines["right"].set_visible(True)
    prev_plot.ax.spines["top"].set_visible(True)

    # horizontal line at 100% prevalence
    prev_plot.ax.axhline(100, color="#666666", ls="--", lw=1.5)

    # adding significance annotations
    chisq_result = chisq(chisq_proc_dataset(data))
    for position, value in enumerate(chisq_result.chisq_result):
        bactprev = plot_data[
            plot_data.cluster == chisq_result.reset_index()["index"][position]
        ].prev_perc
        ypos = max(list(bactprev)) - 1
        if value == "LowN":
            pass
        elif value == "NS":
            pass
        elif value == "P < 0.05":
            prev_plot.ax.text(position, ypos, "*", fontsize=18, ha="center")
            prev_plot.ax.text(position, (ypos + 4), "___", fontsize=18, ha="center")
        elif value == "P < 0.01":
            prev_plot.ax.text(position, ypos, "**", fontsize=18, ha="center")
            prev_plot.ax.text(position, (ypos + 4), "___", fontsize=18, ha="center")
        elif value == "P < 0.001":
            prev_plot.ax.text(position, ypos, "***", fontsize=18, ha="center")
            prev_plot.ax.text(position, (ypos + 4), "___", fontsize=18, ha="center")
        else:
            pass
    return prev_plot


#####
#####           
##########____________________BACTERIOCIN_DIFFERENCES____________________##########
#####
#####

# looking at the differences in prevalence in the pre/post vaccine time periods and between carriage/disease pneumococci 
# plots to visualise the differences, chi square tests for determining which differences are significant 


###
#####_____FIGURE_2B_PRE/POST_PCV10_____#####
###

## specific chi-square function for comparing PCV time periods 
# to be called on either the Icelandic or the Kenyan dataset
# this is fed into the general chisq function 
def chisq_proc_vac(data):
    output = DataFrame(index=unique_clusters(data))
    O1 = []
    O2 = []
    O3 = []
    O4 = []

    for x in unique_clusters(data):
        pre_count = data[data.vaccination == "pre"][
            "{}_status".format(x)
        ].value_counts()
        post_count = data[data.vaccination == "post"][
            "{}_status".format(x)
        ].value_counts()

        if True in pre_count.index:
            O1.append(pre_count[True])
        else:
            O1.append(0)

        if True in post_count.index:
            O2.append(post_count[True])
        else:
            O2.append(0)

        if False in pre_count.index:
            O3.append(pre_count[False])
        else:
            O3.append(0)

        if False in post_count.index:
            O4.append(post_count[False])
        else:
            O4.append(0)

    output["O1"] = O1
    output["O2"] = O2
    output["O3"] = O3
    output["O4"] = O4
    return output

## function which calls the general chisq function and processes the results into a plottable format 
def chisq_results_vac(data):
    # generating an output table to be printed
    # showing the difference in prevalence in the pre/post vac period
    # and showing which differences are statistically signifiacnt using the chisq test

    # setting up output dataframe 
    prev_output = DataFrame(index=unique_clusters(data), columns=["Change_PostPCV10"])

    dataset = list(data.country.unique())
    if len(dataset) > 1:
        return "chisq_results_vac can only be run on one dataset at a time"
    else:
        # setting dataset as a string describing which dataset it has been run on (Iceland or Kenya)
        dataset = dataset[0]

    # increase post vac - positive value
    # postvac - prevac

    for x in prev_output.index:

        pre_count = data[data.vaccination == "pre"][
            "{}_status".format(x)
        ].value_counts()
        post_count = data[data.vaccination == "post"][
            "{}_status".format(x)
        ].value_counts()

        if True in pre_count.index:
            if True in post_count.index:
                prev_output.Change_PostPCV10[x] = (
                    ((post_count[True]) / (sum(post_count))) * 100
                ) - (((pre_count[True]) / (sum(pre_count))) * 100)
            else:
                prev_output.Change_PostPCV10[x] = 0 - (
                    ((pre_count[True]) / (sum(pre_count))) * 100
                )
        else:
            if True in post_count.index:
                prev_output.Change_PostPCV10[x] = (
                    ((post_count[True]) / (sum(post_count))) * 100
                ) - 0
            else:
                prev_output.Change_PostPCV10[x] = 0

    output = pd.concat([prev_output, chisq(chisq_proc_vac(data))], axis=1)
    output["Dataset"] = dataset

    code_col = []
    # populating output with a description of the significance of each difference 
    for x in output.index:
        if output.chisq_result[x] in ["LowN", "NS"]:
            if dataset == "Iceland":
                code_col.append("NS_I")
            else:
                code_col.append("NS_K")
        elif dataset == "Iceland":
            code_col.append("S_I")
        elif dataset == "Kenya":
            code_col.append("S_K")
        else:
            pass
    output["SigCode"] = code_col

    return output

## plotting function 
# using the swarm plot function from seaborn because of the jitter argument - stops the points overlapping
# overlaying three separate scatter plots so that the marker shape can be altered according to the significance of the difference (according to chi sq)
# could use scatterplot, allows marker to be mapped like hue is (for colour), but there is no option for dodge/jitter on a scatter plot, so we can't see all the points
# this is a hacky way to get the figure I want - requires extra data processing, some of the steps are probably redundant, but it works
# final figure is coloured by dataset, non-significant changes in a lighter shade, and different marker indicates degree of significance
def scatter_vac(data):
    I_df = data[data.country == "Iceland"]
    K_df = data[data.country == "Kenya"]
    
    # combine processed dataframes 
    combo = pd.concat([chisq_results_vac(I_df), chisq_results_vac(K_df)])
    combo.reset_index(inplace=True)
    combo.rename(columns={"index": "Bacteriocin"}, inplace=True)

    # function for preparing plot data, called as argument to vac_single_scatter()
    def vac_scatter_prep(sig):
        if sig == "NS":
            temp_data = combo[["Bacteriocin"]].join(
                combo[combo["chisq_result"].isin(["LowN", "NS"])],
                how="left",
                lsuffix="",
                rsuffix="_r",
            )
            temp_data.drop(["Bacteriocin_r"], axis=1, inplace=True)
            return temp_data
        else:
            temp_data = combo[["Bacteriocin"]].join(
                combo[combo["chisq_result"] == "{}".format(sig)],
                how="left",
                lsuffix="",
                rsuffix="_r",
            )
            temp_data.drop(["Bacteriocin_r"], axis=1, inplace=True)
            return temp_data

    # general function for plotting scatters using swarmplot
    def vac_single_scatter(data, size, marker):

        sns.swarmplot(
            data=data,
            x="Bacteriocin",
            y="Change_PostPCV10",
            s=size,
            hue="SigCode",
            marker=marker,
            palette=pal,
            ax=ax1,
            order=x_order,
            linewidth=1,
            edgecolor='black',
        )   
    
    # setting up ax1
    plot, (ax1) = plt.subplots(1, 1, figsize=(8, 4))
    #plt.gcf().subplots_adjust(bottom=0.25)

    # setting up colour palette
    pal = {
        "NS_I": "#85D2FF",
        "NS_K": "#FBBB4B",
        "S_I": sns.color_palette("colorblind")[0],
        "S_K": sns.color_palette("colorblind")[1],
    }
    
    #choosing markers for scatter
    marks = {"LowN": "o", "NS": "o", "P < 0.05": "v", "P < 0.01": "^", "P < 0.001": "d"}

    # manually defining the x axis point order to put the significantly altered ones first
    x_order = ["scy", "sle", "slg", "sce", "scd", "sla", "sls", "slc"]
    for cluster in unique_clusters(data):
        if cluster in x_order:
            pass
        else:
            x_order.append(cluster)

    # plotting scatters
    vac_single_scatter(vac_scatter_prep("NS"), 7, "o")
    vac_single_scatter(vac_scatter_prep("P < 0.001"), 9, "d")
    vac_single_scatter(vac_scatter_prep("P < 0.01"), 9, "^")
    vac_single_scatter(vac_scatter_prep("P < 0.05"), 9, "v")

    # aesthetics
    ax1.grid(axis="both", visible=True)
    ax1.get_legend().set_visible(False)
    ax1.axhline(0, ls="--", color="black")
    ax1.set_ylim(-17, 22)
    ax1.set_xlabel("")
    ax1.set_ylabel("Difference in Prevalence (%)", size=12)
    plt.title("Pre/Post PCV10", fontsize = 12)
    bact_names = [
        "Streptocyclicin",
        "Streptolancidin E",
        "Streptolancidin G",
        "Streptococcin E",
        "Streptococcin D",
        "Streptolancidin A",
        "Streptolassin",
        "Streptolancidin C",
        "Cib",
        "Streptococcin A",
        "Streptococcin B",
        "Streptococcin C",
        "Streptolancidin B",
        "Streptolancidin D",
        "Streptolancidin F",
        "Streptolancidin H",
        "Streptolancidin I",
        "Streptolancidin J",
        "Streptolancidin K",
        "Streptosactin",
    ]
    ax1.set_xticklabels(
        bact_names, rotation=45, size=12, ha="right", rotation_mode="anchor"
    )
    ax1.set_yticklabels(ax1.get_yticks(), size = 12)

    return plot

###
#####_____FIGURE_2C_CARRIAGE/DISEASE_____#####
###

## specific chisq prep function for comparing bacteriocin prevalence in carriage/disease pneumococci 
# written to work on any two categories of disease e.g. carriage vs. IPD, or OM vs. LRTI, two disease arguments this
# uses disease type column: "Carriage", "Invasive", "LRT", "Otitis media"

def chisq_proc_dis(data, disease1, disease2):
    
    dataset = list(data.country.unique())
    if len(dataset) >1:
        return "call chisq_proc_dis on Icelandic or Kenyan data ONLY"
    else:
        dataset = dataset[0]

    data1 = data[data.disease_type == disease1]
    data2 = data[data.disease_type == disease2]

    output = DataFrame(index=unique_clusters(data))

    O1 = []
    O2 = []
    O3 = []
    O4 = []

    for x in unique_clusters(data):

        data1_count = data1["{}_status".format(x)].value_counts()
        data2_count = data2["{}_status".format(x)].value_counts()

        if True in data1_count.index:
            O1.append(data1_count[True])
        else:
            O1.append(0)

        if True in data2_count.index:
            O2.append(data2_count[True])
        else:
            O2.append(0)

        if False in data1_count.index:
            O3.append(data1_count[False])
        else:
            O3.append(0)

        if False in data2_count.index:
            O4.append(data2_count[False])
        else:
            O4.append(0)

    output["O1"] = O1
    output["O2"] = O2
    output["O3"] = O3
    output["O4"] = O4

    return output

## function which calls the general chisq function 
# and wrangles the data into a suitable format for the scatter plot function below 
def chisq_results_dis(data, disease1, disease2):
    # defining a string corresponding to whichever dataset is in use 
    dataset = list(data.country.unique())
    if len(dataset) >1:
        return "call chisq_results_dis on Icelandic or Kenyan data ONLY"
    else:
        dataset = dataset[0]

    data1 = data[data.disease_type == disease1]
    data2 = data[data.disease_type == disease2]

    prev_output = DataFrame(index=unique_clusters(data), columns=["Prev_change"])

    for x in prev_output.index:

        dis1_count = data1["{}_status".format(x)].value_counts()
        dis2_count = data2["{}_status".format(x)].value_counts()

        if True in dis1_count.index:
            if True in dis2_count.index:
                prev_output.Prev_change[x] = (
                    ((dis2_count[True]) / (sum(dis2_count))) * 100
                ) - (((dis1_count[True]) / (sum(dis1_count))) * 100)
            else:
                prev_output.Prev_change[x] = 0 - (
                    ((dis1_count[True]) / (sum(dis1_count))) * 100
                )
        else:
            if True in dis2_count.index:
                prev_output.Prev_change[x] = (
                    ((dis2_count[True]) / (sum(dis2_count))) * 100
                ) - 0
            else:
                prev_output.Prev_change[x] = 0

    output = pd.concat(
        [prev_output, chisq(chisq_proc_dis(data, disease1, disease2))], axis=1
    )
    output["Dataset"] = dataset

    code_col = []
    for x in output.index:
        if output.chisq_result[x] == "LowN":
            if output.Dataset[x] == "Iceland":
                code_col.append("NS_I")
            elif output.Dataset[x] == "Kenya":
                code_col.append("NS_K")
        elif output.chisq_result[x] == "NS":
            if output.Dataset[x] == "Iceland":
                code_col.append("NS_I")
            elif output.Dataset[x] == "Kenya":
                code_col.append("NS_K")
        elif dataset == "Iceland":
            code_col.append("S_I")
        elif dataset == "Kenya":
            code_col.append("S_K")
        else:
            pass
    output["SigCode"] = code_col

    return output

## plotting function, works broadly as scatter_vac (above) 
# option for a 3-figure or 6-figure version using the size argument
# size = "L" returns a figure with additional disease comparisons in the Icelandic dataset (IPD/LRTI etc.)
# size = "" only compares diseases to carriage, used in the paper 

def scatter_dis_panel(data, size):
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    # function for plotting a single panel of the larger scatter plot
    # choose the two subsets and which axis to plot on
    def scatter_dis(data, disease1, disease2, axis):
        
        I_df = data[data.country == "Iceland"]
        K_df = data[data.country == "Kenya"]

        # setting up data to plot - subsets of the appropriate datasets according to arguments to scatter_dis()
        proc_data = "blah"
        if disease1 == "Carriage":
            if disease2 == "Invasive":
                data1 = chisq_results_dis(I_df, "Carriage", "Invasive")
                data2 = chisq_results_dis(K_df, "Carriage", "Invasive")

                proc_data = pd.concat([data1, data2])
            else:
                proc_data = chisq_results_dis(I_df, disease1, disease2)
        else:
            proc_data = chisq_results_dis(I_df, disease1, disease2)

        proc_data.reset_index(inplace=True)
        proc_data.rename(columns={"index": "Bacteriocin"}, inplace=True)

        # general function for preparing data for plotting scatters
        def dis_scatter_prep(sig):
            if sig == "NS":
                temp_data = proc_data[["Bacteriocin"]].join(
                    proc_data[proc_data["chisq_result"].isin(["LowN", "NS"])],
                    how="left",
                    lsuffix="",
                    rsuffix="_r",
                )
                temp_data.drop(["Bacteriocin_r"], axis=1, inplace=True)
                return temp_data

            else:
                temp_data = proc_data[["Bacteriocin"]].join(
                    proc_data[proc_data["chisq_result"] == "{}".format(sig)],
                    how="left",
                    lsuffix="",
                    rsuffix="_r",
                )
                temp_data.drop(["Bacteriocin_r"], axis=1, inplace=True)
                return temp_data

        # general function for plotting scatters
        def dis_single_scatter(data, size, marker):
            sns.swarmplot(
                data=data,
                x="Bacteriocin",
                y="Prev_change",
                s=size,
                hue="SigCode",
                marker=marker,
                palette=pal,
                ax=axis,
                order=x_order,
                linewidth=1,
                edgecolor='black',
            ).set_title(title, size=12)

        # setting up colour palette
        pal = {
            "NS_I": "#85D2FF",
            "NS_K": "#FBBB4B",
            "S_I": sns.color_palette("colorblind")[0],
            "S_K": sns.color_palette("colorblind")[1],
        }

        # generate a title based on the comparison arguments to scatter_dis()

        title = "{} / ".format(disease2) + "{}".format(disease1)

        # manually setting x axis order
        # for now using the order for the car/ipd comparison for all the plots
        x_order = [
            "slc",
            "sls",
            "slf",
            "sca",
            "sce",
            "slg",
            "scd",
            "slb",
            "sle",
            "sla",
            "scy",
            "sld",
        ]
        # adding any bacteriocins not in the manual order in alphabetically afterwards 
        for cluster in unique_clusters(data):
            if cluster in x_order:
                pass
            else:
                x_order.append(cluster)
    
        # calling the scatter plotting functions
        # up to 4 times, once for each possible significance category, using dis_scatter_prep()
        # so that they can be plotted on top of one another for complex aesthetic reasons
        # mostly done within if statements to avoid key errors - not all sig values in all datasets
        dis_single_scatter(dis_scatter_prep("NS"), 7, "o")
        if "P < 0.001" in list(proc_data.chisq_result):
            dis_single_scatter(dis_scatter_prep("P < 0.001"), 9, "d")
        if "P < 0.01" in list(proc_data.chisq_result):
            dis_single_scatter(dis_scatter_prep("P < 0.01"), 9, "v")
        if "P < 0.05" in list(proc_data.chisq_result):
            dis_single_scatter(dis_scatter_prep("P < 0.05"), 9, "^")

    # ____#

    # Setting up axes and dimensions
    if size == "L":
        panel, axs = plt.subplots(2, 3, figsize=(22, 8), sharey=True)
    else:
        panel, axs = plt.subplots(1, 3, figsize=(22, 4), sharey=True)
    axes = axs.flatten()
    plt.subplots_adjust(hspace=0.2, wspace=0.02)

    # plotting scatters
    if size == "L":
        scatter_dis(data, "Carriage", "Invasive", axes[0])
        scatter_dis(data, "Carriage", "Otitis media", axes[1])
        scatter_dis(data, "Carriage", "LRT", axes[2])
        scatter_dis(data, "Invasive", "Otitis media", axes[3])
        scatter_dis(data, "Invasive", "LRT", axes[4])
        scatter_dis(data, "Otitis media", "LRT", axes[5])
    else:
        scatter_dis(data, "Carriage", "Invasive", axes[0])
        scatter_dis(data, "Carriage", "Otitis media", axes[1])
        scatter_dis(data, "Carriage", "LRT", axes[2])

    # x axis labels - full names in order
    bact_names = [
        "Streptolancidin C",
        "Streptolassin",
        "Streptolancidin F",
        "Streptococcin A",
        "Streptococcin E",
        "Streptolancidin G",
        "Streptococcin D",
        "Streptolancidin B",
        "Streptolancidin E",
        "Streptolancidin A",
        "Streptocyclicin",
        "Streptolancidin D",
        "Cib",
        "Streptococcin B",
        "Streptococcin C",
        "Streptolancidin H",
        "Streptolancidin I",
        "Streptolancidin J",
        "Streptolancidin K",
        "Streptosactin",
    ]

    # various aesthetic things
    for x in axes:
        x.grid(axis="both", visible=True)
        x.get_legend().set_visible(False)
        x.set_xticklabels(
            bact_names, rotation=45, size=12, ha="right", rotation_mode="anchor"
        )
        x.axhline(0, ls="--", color="black")
        x.set_ylim(-25, 25)
        x.set_xlabel("")
        x.set_ylabel("")
    axes[0].set_ylabel("Difference in Prevalence (%)", size=12)
    axes[0].set_yticklabels(axes[0].get_yticks(), size = 12)
    if size == "L":
        axes[3].set_ylabel("Difference in Prevalence (%)", size=12)
        for x in [axes[0], axes[1], axes[2]]:
            x.set_xticklabels([])
            x.set_xlabel("")
    else:
        pass
    return panel

#####
#####           
##########____________________BACTERIOCIN/LINEAGE_ASSOCIATIONS____________________##########
#####
#####

###
#####_____GENERAL_FUNCTIONS_____#####
###

## bacteriocin/lineage association function 
# generates a table showing which common lineages a bacteriocin is found in 
# stratified according to a condition (dataset, pre/post vaccination, carriage/disease)
# only shows lineages in which the bacteriocin is detectable at least once 
# if found in >20 lineages, only 20 most common lineages shown, others binned to "other" categories 
# very long and complicated because it works differently depending on the condition 
def bacteriocin_cc_table(data, cluster, condition):
    # data can be full dataset (for location comparison)
    # or Icelandic/Kenyan datasets (for pre/post PCV10 or carriage/disease comparison)

    ##___##
    # setting out the two subsets of the dataset
    # comparing the CC distribution in the two subsets 
    d1 = data
    d2 = data

    # splitting out based on the condition argument 
    if condition == "dataset":
        d1 = data[data.country == "Iceland"]
        d2 = data[data.country == "Kenya"]
        
    elif condition == "vaccine":
        d1 = data[data.vaccination == "pre"]
        d2 = data[data.vaccination == "post"]

    elif condition == "IPD":
        d1 = data[data.carriage_or_disease == "carriage"]
        d2 = data[data.disease_type == "Invasive"]
    
    elif condition == "LRT":
        d1 = data[data.carriage_or_disease == "carriage"]
        d2 = data[data.disease_type == "LRT"]
    
    elif condition == "OM":
        d1 = data[data.carriage_or_disease == "carriage"]
        d2 = data[data.disease_type == "Otitis media"]
        
    else:
        return("No valid condition given as argument")
    ##___##
    

    # the 10 ccs which have the bacteriocin given in the argument most frequently 
    # or all the ccs which have it if it present in less than 10 ccs 
    cc_vcs = data[data["{}_status".format(cluster)] == True].clonal_complex.value_counts()
    ccs = list(data[data["{}_status".format(cluster)] == True].clonal_complex.value_counts().index)
    if len(cc_vcs) > 10:
        ccs = list(data[data["{}_status".format(cluster)] == True].clonal_complex.value_counts().index[:10])
        ccs.append("Other CCs")
        ccs.append("Other Singletons")
    else:
        pass
    
    # start to generate the output data table 
    # this generates the dataframe with the cc column 
    cols = ["CC"]
    cols.append("d1")
    cols.append("d2")
    output = DataFrame(columns = cols)
    
    # populate the lineage column with the relevant ccs 
    output["CC"] = ccs
    
    # need to generate "other" values by iterating through the ccs in the dataframe 
    other_ccs = []
    other_singletons = []
    
    # only doing this if bacteriocin in >10 ccs  
    if "Other CCs" in ccs:
        for cc in list(data[data["{}_status".format(cluster)] == True].clonal_complex.unique()):
            if cc in ccs:
                pass
            elif cc.startswith("S"):
                other_singletons.append(cc)
                # bin to singletons 
            else:
                other_ccs.append(cc)
                # bin to other ccs
    else:
        pass
    
    # function to pool "other" cc totals
    # returns the number of isolates in data that are from ccs in the cc_list
    def other_n(data, cc_list):
        n = 0
        for cc in cc_list:
            if cc in list(data.clonal_complex.unique()):
                n = n + len(data[data.clonal_complex == cc])
            else:
                pass    
        return n 
    
    # returning a string to input directly into the output dataframe 
    # for the other_ccs and other_singletons rows 
    def other_values(data, cluster, cc_list):
        # total number of other CCs 
        total = other_n(data, cc_list)
        if total == 0:
            return "0"
        else:
            # number of other CC isolates with the bacteriocin 
            bact_data = data[data["{}_status".format(cluster)] == True]
            bact_n = other_n(bact_data, cc_list)
            # percentage of other CCs total with the bacteriocin
            bact_perc = (bact_n/total)*100 
            
            if bact_n == 0:
                return "0"
            elif bact_perc == 100:
                value = str(bact_n) + " (" + "100)" 
                return value
            else:
                # manipulate into a string for the output table 
                value = str(bact_n) + " (" + str(round(bact_perc,1)) + ")" 
                return value
    
    # need to fill in the relevant values for each cc
    for row in output.index:
        
        cc = output.CC[row]
        
        if cc == "Other CCs":
            # populating the other CCs field 
            output.d1[row] = other_values(d1, cluster, other_ccs)
            output.d2[row] = other_values(d2, cluster, other_ccs)

        elif cc == "Other Singletons":
            
            # populating the overall count 
            output.d1[row] = other_values(d1, cluster, other_singletons)
            output.d2[row] = other_values(d2, cluster, other_singletons)
            
        else:
            # populating the individual CCs 

            if cc in d1[d1["{}_status".format(cluster)] == True].clonal_complex.value_counts().index:
                total = d1.clonal_complex.value_counts()[cc]

                bact_df = d1[d1.clonal_complex == cc]
                bact_n = len(bact_df[bact_df["{}_status".format(cluster)] == True])
                
                bact_perc = (bact_n/total)*100

                if bact_perc == 100:
                    output.d1[row] = str(bact_n) + " (100)"
                else:
                    output.d1[row] = str(bact_n) + " (" + str(round(bact_perc,1)) + ")"
            else:
                output.d1[row] = "0"

            if cc in d2[d2["{}_status".format(cluster)] == True].clonal_complex.value_counts().index:
                total = d2.clonal_complex.value_counts()[cc]

                bact_df = d2[d2.clonal_complex == cc]
                bact_n = len(bact_df[bact_df["{}_status".format(cluster)] == True])

                bact_perc = (bact_n/total)*100

                if bact_perc == 100:
                    output.d2[row] = str(bact_n) + " (100)"
                else:
                    output.d2[row] = str(bact_n) + " (" + str(round(bact_perc,1)) + ")"
            else:
                output.d2[row] = "0"
            
            # finally adjusting the labels of the lineages
            if cc.startswith("S"):
                pass
            else:
                output.CC[row] = "CC" + cc
            
    return output.set_index("CC")



###
#####_____BETWEEN_DATASETS_____#####
###

##_____TABLES_2_AND_3/SUPPLEMENTARY_TABLE_7____##

## function to generate a lineage association table for each bacteriocin in the given dataset
# comparing lineages in the Kenyan dataset to the Icelandic dataset 
# if an output directory is given as an argument, all the tables will be saved to that directory 
def location_cc_assocations(data, out_dir=None):
    if out_dir:
        print("Saving csvs to {}".format(out_dir))
    else:
        print("No output directory given, printing dataframes for all bacteriocins")
    for cluster in unique_clusters(data):
        output = bacteriocin_cc_table(data, cluster, "dataset")
        if out_dir: 
            directory = out_dir + "{}_lineage_country_comp.csv".format(cluster)
            output.to_csv(directory)
            print("saved csv for {}".format(cluster))
        else:
            print(output)
        

## specific function for breaking down sla distribution in the constituent STs of CC138/176 
# as for above, will only write out to csv if a directory is given 
def cc138_176_table(data, out_dir=None):

    #sts = list(df[df.clonal_complex == "138/176"].mlst_st.value_counts().index)

    cc_df = data[data.clonal_complex == "138/176"]
    
    sts = list(cc_df[cc_df.sla_status == True].mlst_st.value_counts().index)
        
    cols = ["ST", "Iceland_sla_n", "Iceland_n", "Kenya_sla_n", "Kenya_n"]
    output = DataFrame(columns = cols)
    output.ST = sts
    
    I_bact_vcs = cc_df[cc_df.country == "Iceland"][data.sla_status == True].mlst_st.value_counts()
    I_total_vcs = cc_df[cc_df.country == "Iceland"].mlst_st.value_counts()
    K_bact_vcs = cc_df[cc_df.country == "Kenya"][data.sla_status == True].mlst_st.value_counts()
    K_total_vcs = cc_df[cc_df.country == "Kenya"].mlst_st.value_counts()
    
    I_df = data[data.country == "Iceland"]
    K_df = data[data.country == "Kenya"]

    for row in output.index:
        st = output.ST[row]
        if st in I_bact_vcs.index:
            n = I_bact_vcs[st]
            prop = (n/(len(I_df)))*100 
            output.Iceland_sla_n[row] = str(n) + " (" + str(round(prop, 2)) + ")"
        else:
            output.Iceland_sla_n[row] = "0 (0.00)"
        
        if st in I_total_vcs.index:
            n = I_total_vcs[st]
            prop = (n/(len(I_df)))*100 
            output.Iceland_n[row] = str(n) + " (" + str(round(prop, 2)) + ")"
        else:  
            output.Iceland_n[row] = "0 (0.00)"
            
        if st in K_bact_vcs.index:
            n = K_bact_vcs[st]
            prop = (n/(len(K_df)))*100 
            output.Kenya_sla_n[row] = str(n) + " (" + str(round(prop, 2)) + ")"
        else:
            output.Kenya_sla_n[row] = "0 (0.00)"
        
        if st in K_total_vcs.index:
            n = K_total_vcs[st]
            prop = (n/(len(K_df)))*100
            output.Kenya_n[row] = str(n) + " (" + str(round(prop, 2)) + ")"
        else:  
            output.Kenya_n[row] = "0 (0.00)"
    
        output.ST[row] = "ST" + str(st)
    
    if out_dir:
        directory = out_dir + "sla_cc138176_country_comp.csv"
        output.to_csv(directory)
        print("Saved sla/cc138_176 assocation table as csv to " + directory)
    else:
        print("No output directory given, printing sla/cc138_176 association table")
        return output

###
#####_____PRE/POST-PCV10_____#####
###

###
#####_____CARRIAGE/DISEASE_____#####
###

# data are presented together as bar plots (for streptolancidin C and streptocyclicin) and tables (for other sig. diff. bacteriocins)
# bar plots are complex - stacked to show the proportion of isolates of each lineage which had the bacteriocin 
# also grouped within CC but stratified according to pre/post PCV10 time period or carriage/disease 

## data processing for plots: 
# cclist generates a list of ccs which contain the bacteriocin of interest 
# bins any cc which represents < 1% of the dataset to the "other" category 
# cutoff argument can be adjusted to change the % of the cutoff - 1% by default 
def cclist(bacteriocin, data, cutoff=None):
    cc_list = (
        data[data["{}_status".format(bacteriocin)] == True]
    ).clonal_complex.unique()
    output = []
    counts = data.clonal_complex.value_counts(normalize=True)
    # only including ccs which represent more than 1% of the population
    for x in cc_list:
        if cutoff:
            if counts[x] > cutoff:
                output.append(x)
            else:
                pass
        else:
            if counts[x] > 0.01:
                output.append(x)
            else:
                pass

    return output

# cc_ser_axisorder returns a list of ccs in order of overall frequency, which is used in plotting to determine the y axis order 
def cc_ser_axisorder(bacteriocin, data, labellist):
    def bactCC_axisorder(bacteriocin, data):
        CClist = cclist(bacteriocin, data)

        subdata = data[data.clonal_complex.isin(CClist)]
        counts = subdata.clonal_complex.value_counts()
        templist = list(counts.index)
        templist.append("Other")

        return templist

    CClist = bactCC_axisorder(bacteriocin, data)
    yaxisorder = []
    for cc in CClist:
        for cc_ser in labellist:
            if cc_ser.split(" ")[0] == cc:
                yaxisorder.append(cc_ser)
            else:
                pass
    return yaxisorder

## function returning a list of cc names modified to include serotype 
# manually coded in the ones which have >1 dominant serotype 
# plot_data is plot_data within the bact_cc_vac/disease data prep functions (where dom_ser is called)
# data is the overall dataframe, inherited from the bact_cc_vac/disease data prep functions
def dom_ser(data, plot_data):
    
    # manually written in lists for ccs with multiple dominant serotypes
    I_multipleser = ["30", "15"]
    K_multipleser = ["230", "5329", "701", "854", "5902"]
    
    dom_ser_list = []
    for x in plot_data.index:
        tempcc = plot_data.clonal_complex[x]
        tempdf = data[data.clonal_complex == tempcc]
        sercount = tempdf.final_serotype_designation.value_counts()

        if list(data.country.unique())[0] == "Iceland":
            if plot_data.clonal_complex[x] in I_multipleser:
                dom1 = sercount.index[0]
                dom2 = sercount.index[1]
                dom_ser_list.append(tempcc + " (" + dom1 + ", " + dom2 + ")")
            else:
                dom = sercount.index[0]
                dom_ser_list.append(tempcc + " (" + dom + ")")

        elif list(data.country.unique())[0] == "Kenya":
            if plot_data.clonal_complex[x] in K_multipleser:
                dom1 = sercount.index[0]
                dom2 = sercount.index[1]
                dom_ser_list.append(tempcc + " (" + dom1 + ", " + dom2 + ")")
            else:
                dom = sercount.index[0]
                dom_ser_list.append(tempcc + " (" + dom + ")")
    return dom_ser_list

## returns a dict used by all the below processing functions 
# avoiding repeated code
def other_dict(bacteriocin, condition, subset_length, other, other_bact):
    output = {
        "clonal_complex": "Other",
        "CC_ser": "Other",
        "cc_count": len(other),
        "condition": condition,
        "cc_proportion": (len(other) / subset_length) * 100,
        "{}_proportion".format(bacteriocin): (len(other_bact) / subset_length) * 100,
        }
    return output

## processing data for plotting pre/post vac associations 
# complex, and could be more generalised with the other processing functions for the carriage/disease comparisons 
def bact_cc_vac(bacteriocin, data):
    # setting denominators for proportion calculation
    prelength = len(data[data.vaccination == "pre"])
    postlength = len(data[data.vaccination == "post"])
    # generating list of CCs which are present >1% of overall dataset
    CClist = cclist(bacteriocin, data)
    # generate sub-dataframes of only the CCs in cclist both pre and post vaccination
    subdata = data[data.clonal_complex.isin(CClist)]
    pre_subdata = subdata[subdata.vaccination == "pre"]
    post_subdata = subdata[subdata.vaccination == "post"]

    def output_generation(bacteriocin, dataframe):
        count = DataFrame(dataframe.clonal_complex.value_counts())
        count.reset_index(inplace=True)
        count.rename(
            columns={"index": "clonal_complex", "clonal_complex": "cc_count"},
            inplace=True,
        )

        denom = "blah"
        if dataframe is pre_subdata:
            denom = prelength
            count["vaccination"] = len(count) * ["pre"]
        elif dataframe is post_subdata:
            denom = postlength
            count["vaccination"] = len(count) * ["post"]
        else:
            pass

        ccprop = []
        for x in count.index:
            prop = ((count.cc_count[x]) / denom) * 100
            ccprop.append(prop)
        count["cc_proportion"] = ccprop

        bactprop = []
        bactdata = dataframe[dataframe["{}_status".format(bacteriocin)] == True]
        bact_cc_count = bactdata.clonal_complex.value_counts()

        for x in count.clonal_complex:
            if x in bact_cc_count.index:
                prop = ((bact_cc_count[x]) / denom) * 100
                bactprop.append(prop)
            else:
                bactprop.append(0)

        count["{}_proportion".format(bacteriocin)] = bactprop
        return count

    # setting up the output df by merging output_generation results
    outputdfs = [
        output_generation(bacteriocin, pre_subdata),
        output_generation(bacteriocin, post_subdata),
    ]
    plot_data = pd.concat(outputdfs)
    plot_data.rename(columns = {"vaccination":"condition"}, inplace = True)
    plot_data.reset_index(inplace=True)

    # adding the dominant serotype(s) for each CC
    plot_data["CC_ser"] = dom_ser(data, plot_data)

    # generate an "other" row for both pre and post subdataframes
    # binning all serotypes which represent less than 1% of the overall datasets
    # this cutoff is set to very low so that all CCs are included - hacky workaround for the annoying cclist function 
    all_ccs = cclist(bacteriocin, data, 0.0001)

    other_ccs = []
    for cc in all_ccs:
        if cc in CClist:
            pass
        else:
            other_ccs.append(cc)

    otherccs = data[data.clonal_complex.isin(other_ccs)]

    pre_other = otherccs[otherccs.vaccination == "pre"]
    post_other = otherccs[otherccs.vaccination == "post"]

    pre_other_bact = pre_other[pre_other["{}_status".format(bacteriocin)] == True]
    post_other_bact = post_other[post_other["{}_status".format(bacteriocin)] == True]

    other_pre_dict = other_dict(bacteriocin, "pre", prelength, pre_other, pre_other_bact)
    other_post_dict = other_dict(bacteriocin, "post", postlength, post_other, post_other_bact)

    plot_data = plot_data.append(other_pre_dict, ignore_index=True)
    plot_data = plot_data.append(other_post_dict, ignore_index=True)

    return plot_data

## processing for carriage/IPD comparisons
# similar to above vac comparisons 
def bact_cc_inv(bacteriocin, data):
    # setting denominators for proportion calculation
    carlength = len(data[data.disease_type == "Carriage"])
    invlength = len(data[data.disease_type == "Invasive"])

    # generating list of CCs which are present >1% of overall dataset
    CClist = cclist(bacteriocin, data)

    # generate sub-dataframes of only the CCs in cclist both pre and post vaccination
    subdata = data[data.clonal_complex.isin(CClist)]
    car_subdata = subdata[subdata.disease_type == "Carriage"]
    inv_subdata = subdata[subdata.disease_type == "Invasive"]

    def output_generation(bacteriocin, dataframe):

        count = DataFrame(dataframe.clonal_complex.value_counts())
        count.reset_index(inplace=True)
        count.rename(
            columns={"index": "clonal_complex", "clonal_complex": "cc_count"},
            inplace=True,
        )

        denom = "blah"
        if dataframe is car_subdata:
            denom = carlength
            count["car_inv"] = len(count) * ["car"]
        elif dataframe is inv_subdata:
            denom = invlength
            count["car_inv"] = len(count) * ["inv"]
        else:
            pass

        ccprop = []
        for x in count.index:
            prop = ((count.cc_count[x]) / denom) * 100
            ccprop.append(prop)
        count["cc_proportion"] = ccprop

        bactprop = []
        bactdata = dataframe[dataframe["{}_status".format(bacteriocin)] == True]
        bact_cc_count = bactdata.clonal_complex.value_counts()

        for x in count.clonal_complex:
            if x in bact_cc_count.index:
                prop = ((bact_cc_count[x]) / denom) * 100
                bactprop.append(prop)
            else:
                bactprop.append(0)

        count["{}_proportion".format(bacteriocin)] = bactprop
        # count.rename(columns = {"index":"clonal_complex", "clonal_complex":"cc_count"}, inplace = True)
        return count

    outputdfs = [
        output_generation(bacteriocin, car_subdata),
        output_generation(bacteriocin, inv_subdata),
    ]
    plot_data = pd.concat(outputdfs)
    plot_data.reset_index(inplace=True)
    plot_data.rename(columns ={"car_inv":"condition"}, inplace = True)

    plot_data["CC_ser"] = dom_ser(data, plot_data)

    # generate an "other" row for both subdataframes
    all_ccs = cclist(bacteriocin, data, 0.0001)
    other_ccs = []
    for cc in all_ccs:
        if cc in CClist:
            pass
        else:
            other_ccs.append(cc)
    otherccs = data[data.clonal_complex.isin(other_ccs)]
    car_other = otherccs[otherccs.disease_type == "Carriage"]
    inv_other = otherccs[otherccs.disease_type == "Invasive"]

    car_other_bact = car_other[car_other["{}_status".format(bacteriocin)] == True]
    inv_other_bact = inv_other[inv_other["{}_status".format(bacteriocin)] == True]

    other_cardict = other_dict(bacteriocin, "car", carlength, car_other, car_other_bact)
    other_invdict = other_dict(bacteriocin, "inv", invlength, inv_other, inv_other_bact)

    plot_data = plot_data.append(other_cardict, ignore_index=True)
    plot_data = plot_data.append(other_invdict, ignore_index=True)

    return plot_data

## processing for comparing carriage and non IPD disease specifically in Iceland 
# complex again 
def bact_cc_Idisease(data, bacteriocin):
    data = data[data.country == "Iceland"]
    # setting denominators for proportion calculation
    carlength = len(data[data.disease_type == "Carriage"])
    invlength = len(data[data.disease_type == "Invasive"])
    lrtlength = len(data[data.disease_type == "LRT"])
    omlength = len(data[data.disease_type == "Otitis media"])

    # generating list of CCs which are present >1% of overall dataset
    CClist = cclist(bacteriocin, data)

    # generate sub-dataframes of only the CCs in cclist in all disease states
    subdata = data[data.clonal_complex.isin(CClist)]
    car_subdata = subdata[subdata.disease_type == "Carriage"]
    inv_subdata = subdata[subdata.disease_type == "Invasive"]
    lrt_subdata = subdata[subdata.disease_type == "LRT"]
    om_subdata = subdata[subdata.disease_type == "Otitis media"]

    # to be run on each subdataframe, before stitching back into one big dataframe
    def output_generation(bacteriocin, dataframe):
        count = DataFrame(dataframe.clonal_complex.value_counts())
        count.reset_index(inplace=True)
        count.rename(
            columns={"index": "clonal_complex", "clonal_complex": "cc_count"},
            inplace=True,
        )

        denom = "blah"
        if dataframe is car_subdata:
            denom = carlength
            count["car_dis"] = len(count) * ["car"]
        elif dataframe is inv_subdata:
            denom = invlength
            count["car_dis"] = len(count) * ["inv"]
        elif dataframe is lrt_subdata:
            denom = lrtlength
            count["car_dis"] = len(count) * ["lrt"]
        elif dataframe is om_subdata:
            denom = omlength
            count["car_dis"] = len(count) * ["om"]
        else:
            pass

        ccprop = []
        for x in count.index:
            prop = ((count.cc_count[x]) / denom) * 100
            ccprop.append(prop)
        count["cc_proportion"] = ccprop

        bactprop = []
        bactdata = dataframe[dataframe["{}_status".format(bacteriocin)] == True]
        bact_cc_count = bactdata.clonal_complex.value_counts()

        for x in count.clonal_complex:
            if x in bact_cc_count.index:
                prop = ((bact_cc_count[x]) / denom) * 100
                bactprop.append(prop)
            else:
                bactprop.append(0)

        count["{}_proportion".format(bacteriocin)] = bactprop
        return count

    outputdfs = []
    subdfs = [car_subdata, inv_subdata, lrt_subdata, om_subdata]
    for x in subdfs:
        outputdfs.append(output_generation(bacteriocin, x))
    plot_data = pd.concat(outputdfs)
    plot_data.reset_index(inplace=True)
    plot_data.rename(columns = {"car_dis":"condition"}, inplace=True)

    plot_data["CC_ser"] = dom_ser(data, plot_data)

    # generate an "other" row for both subdataframes
    all_ccs = cclist(bacteriocin, data, 0.0001)
    other_ccs = []
    for cc in all_ccs:
        if cc in CClist:
            pass
        else:
            other_ccs.append(cc)
    otherccs = data[data.clonal_complex.isin(other_ccs)]
    car_other = otherccs[otherccs.disease_type == "Carriage"]
    inv_other = otherccs[otherccs.disease_type == "Invasive"]
    lrt_other = otherccs[otherccs.disease_type == "LRT"]
    om_other = otherccs[otherccs.disease_type == "Otitis media"]

    car_other_bact = car_other[car_other["{}_status".format(bacteriocin)] == True]
    inv_other_bact = inv_other[inv_other["{}_status".format(bacteriocin)] == True]
    lrt_other_bact = lrt_other[lrt_other["{}_status".format(bacteriocin)] == True]
    om_other_bact = om_other[om_other["{}_status".format(bacteriocin)] == True]

    other_cardict = other_dict(bacteriocin, "car", carlength, car_other, car_other_bact)
    other_invdict = other_dict(bacteriocin, "inv", invlength, inv_other, inv_other_bact)
    other_lrtdict = other_dict(bacteriocin, "lrt", lrtlength, lrt_other, lrt_other_bact)
    other_omdict = other_dict(bacteriocin, "om", omlength, om_other, om_other_bact)

    plot_data = plot_data.append(other_cardict, ignore_index=True)
    plot_data = plot_data.append(other_invdict, ignore_index=True)
    plot_data = plot_data.append(other_lrtdict, ignore_index=True)
    plot_data = plot_data.append(other_omdict, ignore_index=True)

    return plot_data

## plotting functions 
# called when generating panels below 
# different condition comparisons have different plotting functions 
def bactcc_vacplot(bacteriocin, data, axis):
    # palette is a list of hex values
    vac_pal = ["#A48AA8", "#57267D", "#D3D3D3", "#808080"]
    # generating the y axis labels 
    ylabels = (bact_cc_vac(bacteriocin, data)).CC_ser.unique()

    sns.barplot(
        y="CC_ser",
        x="cc_proportion",
        data=bact_cc_vac(bacteriocin, data),
        hue="condition",
        palette=vac_pal[2:],
        order=cc_ser_axisorder(bacteriocin, data, ylabels),
        ax=axis,
    )
    sns.barplot(
        y="CC_ser",
        x="{}_proportion".format(bacteriocin),
        data=bact_cc_vac(bacteriocin, data),
        hue="condition",
        palette=vac_pal[0:2],
        order=cc_ser_axisorder(bacteriocin, data, ylabels),
        ax=axis,
    )

def bactcc_invplot(bacteriocin, data, axis):
    # generating the y axis labels   
    ylabels = (bact_cc_inv(bacteriocin, data)).CC_ser.unique()

    sns.barplot(
        y="CC_ser",
        x="cc_proportion",
        data=bact_cc_inv(bacteriocin, data),
        hue="condition",
        palette=["#DCDCDC", "#808080"],
        ax=axis,
        order=cc_ser_axisorder(bacteriocin, data, ylabels),
    )

    sns.barplot(
        y="CC_ser",
        x="{}_proportion".format(bacteriocin),
        data=bact_cc_inv(bacteriocin, data),
        hue="condition",
        palette=[disease_pal[0], disease_pal[1]],
        ax=axis,
        order=cc_ser_axisorder(bacteriocin, data, ylabels),
    )

# this one is more complex as it needs to check which disease states actually show sig differences 
def bactcc_Idisplot(bacteriocin, data, axis):
    # generating a dict of lists describing which subsets of the data each bacteriocin is associated to
    # bacteriocins are associated with different types of disease over carriage in the Icelandic dataset
    subset = ["car"]
    if (
        chisq_results_dis(data[data.country == "Iceland"], "Carriage", "Invasive").SigCode[bacteriocin]
        == "S_I"
    ):
        subset.append("inv")
    if chisq_results_dis(data[data.country == "Iceland"], "Carriage", "LRT").SigCode[bacteriocin] == "S_I":
        subset.append("lrt")
    if (
        chisq_results_dis(data[data.country == "Iceland"], "Carriage", "Otitis media").SigCode[bacteriocin]
        == "S_I"
    ):
        subset.append("om")

    # generating y axis labels 
    ylabels = (bact_cc_Idisease(data, bacteriocin)).CC_ser.unique()

    if subset == ["car", "inv", "lrt", "om"]:
        sns.barplot(
            y="CC_ser",
            x="cc_proportion",
            data=bact_cc_Idisease(bacteriocin),
            hue="condition",
            palette=["#DCDCDC", "#808080", "#C0C0C0", "#A9A9A9"],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )
        sns.barplot(
            y="CC_ser",
            x="{}_proportion".format(bacteriocin),
            data=bact_cc_Idisease(bacteriocin),
            hue="condition",
            palette=[disease_pal[0], disease_pal[1], disease_pal[2], disease_pal[3]],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

    elif subset == ["car", "lrt", "om"]:
        plotdata = bact_cc_Idisease(data, bacteriocin)[
            bact_cc_Idisease(data, bacteriocin).condition.isin(subset)
        ]

        sns.barplot(
            y="CC_ser",
            x="cc_proportion",
            data=plotdata,
            hue="condition",
            palette=["#DCDCDC", "#C0C0C0", "#A9A9A9"],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

        sns.barplot(
            y="CC_ser",
            x="{}_proportion".format(bacteriocin),
            data=plotdata,
            hue="condition",
            palette=[disease_pal[0], disease_pal[2], disease_pal[3]],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

    elif subset == ["car", "inv", "om"]:
        plotdata = bact_cc_Idisease(bacteriocin)[
            bact_cc_Idisease(bacteriocin).condition.isin(subset)
        ]

        sns.barplot(
            y="CC_ser",
            x="cc_proportion",
            data=plotdata,
            hue="condition",
            palette=["#DCDCDC", "#C0C0C0", "#696969"],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

        sns.barplot(
            y="CC_ser",
            x="{}_proportion".format(bacteriocin),
            data=plotdata,
            hue="condition",
            palette=[disease_pal[0], disease_pal[1], disease_pal[3]],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

    elif subset == ["car", "om"]:
        plotdata = bact_cc_Idisease(bacteriocin)[
            bact_cc_Idisease(bacteriocin).condition.isin(subset)
        ]

        sns.barplot(
            y="CC_ser",
            x="cc_proportion",
            data=plotdata,
            hue="condition",
            palette=["#DCDCDC", "#696969"],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )

        sns.barplot(
            y="CC_ser",
            x="{}_proportion".format(bacteriocin),
            data=plotdata,
            hue="condition",
            palette=[disease_pal[0], disease_pal[3]],
            ax=axis,
            order=cc_ser_axisorder(bacteriocin, data, ylabels),
        )


##_____FIGURE_3A_STREPTOLANCIDIN_C____##

## plotting the distribution of slc in genetic lineages of both datasets 
# more common pre-PCV10, and in various disease states relative to carriage, in both Icelandic and Kenyan datasets 
def slc_panel(data):
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    plot = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=5,
        figure=plot,
        width_ratios=[1, 1, 0.75, 1, 1],
        height_ratios=[1],
        wspace=0,
        hspace=0.5,
    )
    ax1 = plot.add_subplot(gs[0, 0])
    ax2 = plot.add_subplot(gs[0, 1])
    ax3 = plot.add_subplot(gs[0, 3])
    ax4 = plot.add_subplot(gs[0, 4])

    axes = [ax1, ax2, ax3, ax4]

    bactcc_vacplot("slc", data[data.country == "Iceland"], ax1)
    bactcc_Idisplot("slc", data[data.country == "Iceland"], ax2)

    bactcc_vacplot("slc", data[data.country == "Kenya"], ax3)
    bactcc_invplot("slc", data[data.country == "Kenya"], ax4)

    for x in axes:
        x.legend_.remove()
        x.set_xlabel("")
        x.set_yticklabels(x.get_yticklabels(), fontsize=12)
        x.set_ylabel("")
        x.xaxis.set_major_locator(MultipleLocator(10))
        x.grid(True)
    
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xticklabels([0, 0, 5, 10, 15])

    for x in [ax1, ax3]:
        x.set_xlabel("% of Pre/Post PCV10\nIsolates", fontsize=12)
    for x in [ax2, ax4]:
        x.set_xlabel("% of Carriage/Disease\nIsolates", fontsize=12)
        x.set_yticklabels([])

    ax1.set_ylabel("Clonal Complex (Predominant Serotype)", fontsize=12)

    ax2.set_title("Iceland\nCarriage/Disease", fontsize=12)
    ax1.set_title("Iceland\nPre/Post PCV10", fontsize=12)
    ax4.set_title("Kenya\nCarriage/Disease", fontsize=12)
    ax3.set_title("Kenya\nPre/Post PCV10", fontsize=12)

    return plot

##_____FIGURE_3B_STREPTOCYCLICIN____##

def scy_panel(data):
    plot = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=5,
        figure=plot,
        width_ratios=[1, 1, 0.75, 1, 1],
        height_ratios=[1],
        wspace=0,
        hspace=0.5,
    )

    ax1 = plot.add_subplot(gs[0, 0])
    ax2 = plot.add_subplot(gs[0, 1])
    ax3 = plot.add_subplot(gs[0, 3])
    ax4 = plot.add_subplot(gs[0, 4])

    axes = [ax1, ax2, ax3, ax4]

    bactcc_vacplot("scy", data[data.country == "Iceland"], ax1)
    bactcc_Idisplot("scy", data[data.country == "Iceland"], ax2)

    bactcc_vacplot("scy", data[data.country == "Kenya"], ax3)
    bactcc_invplot("scy", data[data.country == "Kenya"], ax4)

    for x in axes:
        x.legend_.remove()
        x.set_xlabel("")
        x.set_ylabel("")
        x.set_yticklabels(x.get_yticklabels(), fontsize=12)
        x.grid(True)

    ax1.set_ylabel("Clonal Complex (Predominant Serotype)", fontsize=12)
    ax3.set_xticklabels([0, 5, 10])

    for x in [ax1, ax3]:
        x.set_xlabel("% of Pre/Post PCV10\nIsolates", fontsize=12)
    for x in [ax2, ax4]:
        x.set_xlabel("% of Carriage/Disease\nIsolates", fontsize=12)
        x.set_yticklabels([])

    ax2.set_title("Iceland\nCarriage/Disease", fontsize=12)
    ax1.set_title("Iceland\nPre/Post PCV10", fontsize=12)
    ax4.set_title("Kenya\nCarriage/Disease", fontsize=12)
    ax3.set_title("Kenya\nPre/Post PCV10", fontsize=12)

    return plot


##_____SUPPLEMENTARY_TABLES_8/9_____##

## generates a table for all relevant vaccine/disease categories for a single bacteriocin 
# works on a single dataset at a time because the Icelandic and Kenyan disease processes differ 
# saves out as csv if argument given 
# manually join and format whichever tables are required in Excel/Word 
def bacteriocin_cc_table_combo(data, cluster):
    
    dataset = "blah"
    # different for Iceland and Kenya - first figure out which dataset it is 
    if len(data.country.unique()) > 1:
        print("Only run on one dataset at a time")
        return 
    elif "Iceland" in data.country.unique():
        output = pd.concat([bacteriocin_cc_table(data, cluster, "vaccine").rename(columns={"d1": "Pre-PCV",  "d2": "Post-PCV"}),
                            bacteriocin_cc_table(data, cluster, "IPD").rename(columns={"d1": "Carriage", "d2": "IPD"}), 
                           bacteriocin_cc_table(data, cluster, "LRT").rename(columns={"d1": "Carriage",  "d2": "LRTI"}).drop(['Carriage'], axis=1), 
                           bacteriocin_cc_table(data, cluster, "OM").rename(columns={"d1": "Carriage",  "d2": "OM"}).drop(['Carriage'], axis=1)], axis = 1)

    elif "Kenya" in data.country.unique():
        output = pd.concat([bacteriocin_cc_table(data, cluster, "vaccine").rename(columns={"d1": "Pre-PCV",  "d2": "Post-PCV"}), 
                            bacteriocin_cc_table(data, cluster, "IPD").rename(columns={"d1": "Carriage", "d2": "IPD"})], axis = 1)

        
    return output 

## function to generate and save combo tables for each cluster in each dataset 
# each table includes lineages split out by pre/post PCV10 and carriage/disease
# generates a full table for every bacteriocin observed in each dataset, even if no significantly different prevalences observed
# choose the ones of interest by hand afterwards 
def bacteriocin_cc_table_combo_outputs(data, country, out_dir):
    subdata = data[data.country == country]
    for cluster in unique_clusters(subdata):
        bacteriocin_cc_table_combo(subdata, cluster).to_csv((out_dir) + "_{}_".format(country) + "{}.csv".format(cluster))
        print("Saved csv for {}".format(cluster) + " in {}".format(country))


    
###
#####_____BACTERIOCIN/SEROTYPE_ASSOCIATIONS_____#####
###

# as above for bacteriocin/CC associations, want to see how the bacteriocins are distributed in the serotypes of the datasets 
# only doing this for the bacteriocins which are more common in a disease process than in carriage 

def bacteriocin_serotype_table(data, cluster):
    # data should be I_df or K_df 
    
    # the 10 serotypes which have the bacteriocin given in the argument most frequently 
    # or all the serotypess which have it if it present in less than 10
    serotype_vcs = data[data["{}_status".format(cluster)] == True].final_serotype_designation.value_counts()
    serotypes = list(serotype_vcs.index)
    if len(serotype_vcs) > 10:
        serotypes = list(serotype_vcs.index[:10])
        serotypes.append("Other Serotypes")
    else:
        pass
    
    # start to generate the output data table 
    # this generates the dataframe with the serotype column 
    cols = ["Serotype"]
    cols.append("n")
    output = DataFrame(columns = cols)
    
    # populate the serotype column with the relevant serotypes
    output["Serotype"] = serotypes
    
    # need to generate "other" values by iterating through the ccs in the dataframe 
    other_serotypes = []
    
    # only doing this if bacteriocin in >10 ccs  
    if "Other Serotypes" in serotypes:
        for serotype in list(data[data["{}_status".format(cluster)] == True].final_serotype_designation.unique()):
            if serotype in serotypes:
                pass
            else:
                other_serotypes.append(serotype)
                # bin to other serotypes
    else:
        pass
    
    # function to pool "other" serotype totals
    # returns the number of isolates in data that are from serotypes in the serotype_list
    def other_n(data, serotype_list):
        n = 0
        for serotype in serotype_list:
            if serotype in list(data.final_serotype_designation.unique()):
                n = n + len(data[data.final_serotype_designation == serotype])
            else:
                pass    
        return n 

    # returning a string to input directly into the output dataframe 
    # for the other_serotypes row 
    def other_values(data, cluster, serotype_list):
        # total number of other CCs 
        total = other_n(data, serotype_list)
        if total == 0:
            return "0"
        else:
            # number of other serotype isolates with the bacteriocin 
            bact_data = data[data["{}_status".format(cluster)] == True]
            bact_n = other_n(bact_data, serotype_list)
            # percentage of other serotypes total with the bacteriocin
            bact_perc = (bact_n/total)*100 
            
            if bact_perc == 100:
                value = str(bact_n) + " (100)" 
                return value
            else:
                value = str(bact_n) + " (" + str(round(bact_perc,1)) + ")" 
                return value
    
    # fill in the relevant values for each serotype
    for row in output.index:
        
        serotype = output.Serotype[row]
        
        if serotype == "Other Serotypes":
            # populating the other serotypes field 
            output.n[row] = other_values(data, cluster, other_serotypes)
          
        else:
            # populating the individual serotypes 
            # total number of examples of the serotype in the dataset 
            total = data.final_serotype_designation.value_counts()[serotype]
            # total number of the serotype which have the bacteriocin 
            bact_df = data[data.final_serotype_designation == serotype]
            bact_n = len(bact_df[bact_df["{}_status".format(cluster)] == True])
            
            # expressed as a %
            bact_perc = (bact_n/total)*100

            if bact_n == 0:
                output.n[row] = "0"
            elif bact_perc == 100:
                output.n[row] = str(bact_n) + " (100)"
            else:
                output.n[row] = str(bact_n) + " (" + str(round(bact_perc,1)) + ")"
    
    output.rename(columns = {"n":"n (%)"}, inplace = True)

    return output.set_index("Serotype")

# function to save csvs for all the bacteriocins showing their serotype associations in the given dataset 
def bacteriocin_serotype_table_outputs(data, output_dir):

    for cluster in unique_clusters(data):

        output = bacteriocin_serotype_table(data, cluster)

        output.to_csv((output_dir + "_" + cluster + ".csv"))



#####
#####           
##########____________________BACTERIOCIN_REPERTOIRES____________________##########
#####
#####

###
#####_____CLUSTERS_PER_GENOME_____#####
###

# how many bacteriocin clusters are present in each genome?
# using the cluster_count field in the processed dataframe to plot the frequency of each number of bacteriocins 

##_____FIGURE_4A_____##

# Uses the seaborn count plot and takes data straight from the dataframe
# This is technically a histogram - could adjust the bar spacing to remove gaps between categories

def cluster_count_bar(data):

    sns.set_style("ticks", {"axes.grid": True})
    sns.set_context("paper")

    count_bar = sns.catplot(
        x="cluster_count",
        # row = "country",
        hue="country",
        kind="count",
        data=data,
        dodge=True,
        legend=False,
        height=4.5,
        aspect=1.5,
    )

    axes = count_bar.axes.flatten()
    axes[0].set_xlabel("Clusters per Genome", fontsize=12)
    axes[0].set_ylabel("Frequency", fontsize=12)
    axes[0].set_xticklabels(axes[0].get_xticklabels(), fontsize=12)
    axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=12)
    axes[0].spines["right"].set_visible(True)
    axes[0].spines["top"].set_visible(True)

    axes[0].legend(loc="upper right", title="Dataset", fontsize=12, title_fontsize=12)

    # plt.subplots_adjust(top=0.85)
    # count_bar.fig.suptitle('Bacteriocin Clusters Detected per Genome')

    return count_bar

###
#####_____BACTERIOCIN_REPERTOIRE_____#####
###

## Repertoires - general heatmap functions 

def general_heatmap(data, title, procdf):
    # data is an array of arrays used for the heatmap itself
    # title is formatted in the processing function, a string describing the heatmap
    # procdf is the processed dataframe generated in the processing function, used for annotating the heatmap

    # set up axis
    height = (len(procdf)) / 5
    plot, (ax1) = plt.subplots(1, figsize=(8, height))
    plt.subplots_adjust(top=0.88, hspace=0.5, wspace=0.5) 
    
    # defining a color map to just give black/grey for presence/absence of each bacteriocin 
    mono_cmap = ["#FFFFFF"] + 20*["#708090"]
    #["#2F4F4F"]
    
    # setting up x axis labels
    bact_names = [
        "Cib",
        "Streptococcin A",
        "Streptococcin B",
        "Streptococcin C",
        "Streptococcin D",
        "Streptococcin E",
        "Streptocyclicin",
        "Streptolancidin A",
        "Streptolancidin B",
        "Streptolancidin C",
        "Streptolancidin D",
        "Streptolancidin E",
        "Streptolancidin F",
        "Streptolancidin G",
        "Streptolancidin H",
        "Streptolancidin I",
        "Streptolancidin J",
        "Streptolancidin K",
        "Streptolassin",
        "Streptosactin",
    ]

    # plotting the heatmap
    plot = sns.heatmap(
        data,
        vmin=0,
        vmax=20,
        cmap=mono_cmap,
        cbar=False,
        linewidths=0.3,
        yticklabels=False,
        ax=ax1,
        linecolor = "black", 
        linewidth = 0.3
    )

    # add title and x tick labels
    ax1.set_title("{}".format(title), fontsize=12, loc="left")
    ax1.set_xticklabels(
        bact_names, fontsize=12, rotation=30, ha="right", rotation_mode="anchor"
    )

    # plotting on the frequencies and proportions
    for position, value in enumerate(list(procdf["count"])):
        prop = procdf.proportion[position]
        freq = str(value) + " (" + str(prop) + ")"
        ax1.text(20.5, (position + 0.7), freq, fontsize=12)
    ax1.text(20.5, -0.5, "n (%)", fontsize=12)

    # adding an * next to any repertoire which is unique to the dataset in question
    # only used in Iceland/Kenya comparison
    if "unique" in procdf.columns:
        for position, value in enumerate(list(procdf.unique)):
            astx = procdf.unique[position]
            ax1.text(20, (position + 1.5), astx, fontsize=20)
    else:
        pass

    return plot



##_____FIGURE_4B_____##

# Heatmaps showing the top 10 most common repertoires in each dataset 
# Previous versions of this include option to plot subsets (carriage/disease or pre/postPCV), removed here because not in use 

# processing data for plotting the top 10 repertoires from each dataset
def bact_rep_top10proc(data, country):

    def shared_rep():
        shared_reps = []
        for rep in list(data[data.country == "Iceland"].bacteriocin_repertoire.unique()):
            if rep in list(data[data.country == "Kenya"].bacteriocin_repertoire.unique()):
                shared_reps.append(rep)
            else:
                pass

        unique_col = []

        for x in top10.index:
            if top10.bacteriocin_repertoire[x] in shared_reps:
                unique_col.append("")
            else:
                unique_col.append("*")
        top10["unique"] = unique_col

    # data prep
    top10 = DataFrame(data[data.country == country].bacteriocin_repertoire.value_counts()[0:10])
    top10.reset_index(inplace=True)
    top10.rename(
        columns={"bacteriocin_repertoire": "count", "index": "bacteriocin_repertoire"},
        inplace=True,
    )

    proplist = []
    for x in top10.index:
        prop = ((top10["count"][x]) / (len(data))) * 100
        prop = round(prop, 1)
        proplist.append(str(prop))
    top10["proportion"] = proplist

    # call the shared_rep function only if looking at overall top10 in each dataset
    shared_rep()

    # generating the array data used to plot the heatmaps
    heatmap_array = []

    for rep in top10.index:
        bact_rep = top10.bacteriocin_repertoire[rep]
        bact_rep_list = bact_rep.split("-")

        single_array = []

        for position, bacteriocin in enumerate(bact_rep_list):
            if bacteriocin == "/":
                single_array.append(0)
            else:
                single_array.append(position + 1)

        heatmap_array.append(single_array)

    title = "{}, Top 10 Repertoires".format(country)

    return heatmap_array, title, top10

# plotting overall top 10 repertoires in both datasets
def top10rep_hm(data):
    plot1 = general_heatmap(*bact_rep_top10proc(data, "Iceland"))
    plot1.set_xticklabels([])
    plot2 = general_heatmap(*bact_rep_top10proc(data, "Kenya"))

    return plot1.get_figure(), plot2.get_figure()


###___MIXED_REPERTOIRES_WITHIN_CC___###

# looking a mixed repertoires within CCs and STs
# first some dscriptive functions returning text lists of CCs with mixed repertoires 

def mixed_rep_cc(data):

    # run this one on one dataset at a time
    mixed_ccs = []

    for CC in list(data.clonal_complex.value_counts().index):
        temp = data[data.clonal_complex == CC]
        if (
            len(temp.groupby(["clonal_complex"]).bacteriocin_repertoire.value_counts())
            > 1
        ):
            mixed_ccs.append(CC)
        else:
            pass

    # this list is sorted in order of most to least frequent
    return mixed_ccs

def cc_repertoires(data):

    print("Iceland\n")
    print("Number of clonal complexes with mixed bacteriocin repertoires:")
    print(len(mixed_rep_cc(data[data.country == "Iceland"])))

    print("\nKenya\n")
    print("Number of clonal complexes with mixed bacteriocin repertoires:")
    print(len(mixed_rep_cc(data[data.country == "Kenya"])))

# function used to examine repertoires of STs within CCs with mixed repertoires 
def mixedcc_stcheck(data):

    result_dict = {}

    # dict with key of ccs with sts with mixing
    # list of all sts with mixed bacteriocin repertoires from that cc

    for cc in mixed_rep_cc(data):

        st_list = []
        tempdf = data[data.clonal_complex == cc]

        for st in tempdf.mlst_st.unique():
            st_subset = (
                tempdf[tempdf.mlst_st == st]
                .groupby(["mlst_st"])
                .bacteriocin_repertoire.value_counts()
            )
            st_subset_df = DataFrame(st_subset)
            st_subset_df.rename(
                columns={"bacteriocin_repertoire": "count"}, inplace=True
            )
            st_subset_df.reset_index(inplace=True)
            # print(st_subset_df)

            if len(st_subset) > 1:
                st_list.append(st)

        if len(st_list) > 0:
            result_dict[cc] = st_list

    return result_dict

# returns descriptive text of how many STs in each dataset have mixed bacteriocin repertoires 
def mixed_sts():
    
    I_total = 0
    for cc in mixedcc_stcheck("Iceland").keys():
        I_total = I_total + len(mixedcc_stcheck("Iceland")[cc])
    
    K_total = 0
    for cc in mixedcc_stcheck("Kenya").keys():
        K_total = K_total + len(mixedcc_stcheck("Kenya")[cc])

    print("Number of STs with multiple bacteriocin repertoires in the Icelandic dataset:")
    print(I_total)
    print("Number of STs with multiple bacteriocin repertoires in the Kenyan dataset:")
    print(K_total)

## Functions to generate a table describing the CCs with mixed repertoires 

def lineage_rep_table(data):

    bact_names = {
        "cib": "Cib",
        "sca": "Streptococcin A",
        "scb": "Streptococcin B",
        "scc": "Streptococcin C",
        "scd": "Streptococcin D",
        "sce": "Streptococcin E",
        "scy": "Streptocyclicin",
        "sla": "Streptolancidin A",
        "slb": "Streptolancidin B",
        "slc": "Streptolancidin C",
        "sld": "Streptolancidin D",
        "sle": "Streptolancidin E",
        "slf": "Streptolancidin F",
        "slg": "Streptolancidin G",
        "slh": "Streptolancidin H",
        "sli": "Streptolancidin I",
        "slj": "Streptolancidin J",
        "slk": "Streptolancidin K",
        "sls": "Streptolassin",
        "ssa": "Streptosactin",
    }

    ccdf_list = []

    # iterate through the CCs which have multiple bacteriocin repertoires
    for cc in mixed_rep_cc(data):
        ccdf = DataFrame(
            columns=[
                "CC",
                "Variable bacteriocins (CC)",
                "Mixed STs",
                "Variable bacteriocins (ST)",
            ]
        )

        cc_list = []
        var_bact_cc = []
        st_list = []
        var_bact_st = []

        # generate a list of 20 empty lists
        rep_list = []
        for i in range(0, 20):
            rep_list.append([])

        mixed_list = []
        mixed_string = ""

        # populate rep_list
        for rep in list(
            data[data.clonal_complex == cc].bacteriocin_repertoire.unique()
        ):
            for position, value in enumerate(list(rep.split("-"))):
                if value in rep_list[position]:
                    pass
                else:
                    rep_list[position].append(value)

        # use rep_list to populate mixed_list
        for x in rep_list:
            if len(x) == 1:
                pass
            else:
                for value in x:
                    if value == "/":
                        pass
                    else:
                        mixed_list.append(value)
        #print(rep_list)
        #print(mixed_list)
        # use mixed_list to generate string of the variable bacteriocins in the CC
        for x in mixed_list:
            if mixed_string == "":
                mixed_string = mixed_string + bact_names[x]
            else:
                mixed_string = mixed_string + ", " + bact_names[x]

        ###
        # same again for bacteriocins with ST
        if cc in mixedcc_stcheck(data).keys():
            for st in mixedcc_stcheck(data)[cc]:
                cc_list.append(cc)
                var_bact_cc.append(mixed_string)
                st_list.append(st)

                mixed_st_list = []
                mixed_st_string = ""

                # generate a list of 20 empty lists
                rep_st_list = []
                for i in range(0, 20):
                    rep_st_list.append([])

                # populate rep_st_list
                for rep in list(
                    data[data.mlst_st == st].bacteriocin_repertoire.unique()
                ):
                    for position, value in enumerate(list(rep.split("-"))):
                        if value in rep_st_list[position]:
                            pass
                        else:
                            rep_st_list[position].append(value)

                # use rep_st_list to populate mixed_list
                for x in rep_st_list:
                    if len(x) == 1:
                        pass
                    else:
                        for value in x:
                            if value == "/":
                                pass
                            else:
                                mixed_st_list.append(value)

                # use mixed_st_list to generate string of the variable bacteriocins in the CC
                for x in mixed_st_list:
                    if mixed_st_string == "":
                        mixed_st_string = mixed_st_string + bact_names[x]
                    else:
                        mixed_st_string = mixed_st_string + ", " + bact_names[x]

                var_bact_st.append(mixed_st_string)
        else:
            cc_list.append(cc)
            var_bact_cc.append(mixed_string)
            st_list.append("None")
            var_bact_st = "NA"

        ccdf["CC"] = cc_list
        ccdf["Variable bacteriocins (CC)"] = var_bact_cc
        ccdf["Mixed STs"] = st_list
        ccdf["Variable bacteriocins (ST)"] = var_bact_st

        ccdf_list.append(ccdf)

    output = pd.concat(ccdf_list)
    return output.reset_index().drop(["index"], axis=1)



##_____FIGURE_4C_____##
##_____FIGURE_4D_____##


## Processing functions 

# call on a single country at a time 
# e.f. df[df.country == "Iceland"]
def bact_rep_ccproc(data, cc):

    # data prep
    cc_subset = (
        data[data.clonal_complex == cc]
        .groupby(["clonal_complex"])
        .bacteriocin_repertoire.value_counts()
    )
    cc_subset_df = DataFrame(cc_subset)
    cc_subset_df.rename(columns={"bacteriocin_repertoire": "count"}, inplace=True)
    cc_subset_df.reset_index(inplace=True)

    proplist = []
    for x in cc_subset_df.index:
        prop = ((cc_subset_df["count"][x]) / (sum(cc_subset_df["count"]))) * 100
        prop = round(prop, 1)
        proplist.append(prop)

    cc_subset_df["proportion"] = proplist

    # print(proplist)
    # print(cc_subset_df)

    heatmap_array = []

    for x in cc_subset_df.index:
        bact_rep = cc_subset_df.bacteriocin_repertoire[x]
        bact_rep_list = bact_rep.split("-")

        single_array = []

        for position, bacteriocin in enumerate(bact_rep_list):
            if bacteriocin == "/":
                single_array.append(0)
            else:
                single_array.append(position + 1)

        heatmap_array.append(single_array)

    title = "{}, ".format(data.country.unique()[0]) + "CC{}".format(cc)

    return heatmap_array, title, cc_subset_df

# call on a single country at a time 
# e.f. df[df.country == "Iceland"]
def bact_rep_stproc(data, st):

    # data prep
    st_subset = (
        data[data.mlst_st == st]
        .groupby(["mlst_st"])
        .bacteriocin_repertoire.value_counts()
    )
    st_subset_df = DataFrame(st_subset)
    st_subset_df.rename(columns={"bacteriocin_repertoire": "count"}, inplace=True)
    st_subset_df.reset_index(inplace=True)

    proplist = []
    for x in st_subset_df.index:
        prop = ((st_subset_df["count"][x]) / (sum(st_subset_df["count"]))) * 100
        prop = round(prop, 1)
        proplist.append(prop)

    st_subset_df["proportion"] = proplist

    # print(proplist)
    # print(cc_subset_df)

    heatmap_array = []

    for x in st_subset_df.index:
        bact_rep = st_subset_df.bacteriocin_repertoire[x]
        bact_rep_list = bact_rep.split("-")

        single_array = []

        for position, bacteriocin in enumerate(bact_rep_list):
            if bacteriocin == "/":
                single_array.append(0)
            else:
                single_array.append(position + 1)

        heatmap_array.append(single_array)

    title = "{}, ".format(data.country.unique()[0]) + "ST{}".format(st)

    return heatmap_array, title, st_subset_df


## Plotting functions 

# Each return a set of heatmaps as separate plots to be stitched together  

def CC439_hm(data):
    # only care about CC439 in the Icelandic dataset 
    # this is here just in case the function is called on the whole dataset 
    sub_data = data[data.country == "Iceland"]
    plot_dict = {
        "cc439": general_heatmap(*bact_rep_ccproc(sub_data, "439")),
        "st311": general_heatmap(*bact_rep_stproc(sub_data, 311)),
        "st507": general_heatmap(*bact_rep_stproc(sub_data, 507)),
        "st190": general_heatmap(*bact_rep_stproc(sub_data, 190)),
        "st442": general_heatmap(*bact_rep_stproc(sub_data, 442)),
    }

    for x in plot_dict.keys():
        if x == "cc439":
            pass
        elif x == "st442":
            pass
        else:
            plot_dict[x].set_xticklabels([])

        plot_dict[x] = plot_dict[x].get_figure()

    return plot_dict


def CC5902_hm(data):
    # only care about CC5902 in the Kenyan dataset 
    # this is here just in case the function is called on the whole dataset 
    sub_data = data[data.country == "Kenya"]
    plot_dict = {
        "cc5902": general_heatmap(*bact_rep_ccproc(sub_data, "5902")),
        "st5370": general_heatmap(*bact_rep_stproc(sub_data, 5370)),
        "st840": general_heatmap(*bact_rep_stproc(sub_data, 840)),
        "st2052": general_heatmap(*bact_rep_stproc(sub_data, 2052)),
        "st5902": general_heatmap(*bact_rep_stproc(sub_data, 5902)),
        "st15056": general_heatmap(*bact_rep_stproc(sub_data, 15056)),
    }

    for x in plot_dict.keys():
        if x == "cc5902":
            pass
        elif x == "st15056":
            pass
        else:
            plot_dict[x].set_xticklabels([])

        plot_dict[x] = plot_dict[x].get_figure()

    return plot_dict



#####
#####           
##########____________________SEROTYPE_PLOT____________________##########
#####
#####

# This plot is not included in the paper, but is included in the thesis chapter 

##_____PROCESSING_____##

def serotype_process(data):  # data is either I_df or K_df, output should be directly plottable

    # generating the grouped data
    # big_df has all the serotypes, before being binned into other or counted for stacked bars

    xs = data.final_serotype_designation.unique()
    ys = data.vaccination.unique()
    zs = data.disease_type.unique()

    ind = pd.MultiIndex.from_product(
        [xs, ys, zs],
        names=["final_serotype_designation", "vaccination", "disease_type"],
    )
    a = np.zeros(len(ind))

    big_df = pd.DataFrame(a, index=ind)
    big_df.reset_index(inplace=True)

    ser = big_df.final_serotype_designation
    ser = ser.str.cat(big_df.vaccination, sep="_")

    big_df["sercat"] = ser
    big_df.drop(0, axis=1, inplace=True)
    big_df.set_index(
        ["final_serotype_designation", "vaccination", "disease_type"], inplace=True
    )
    big_df = big_df.reindex(
        ["Carriage", "Invasive", "LRT", "Otitis media"], level="disease_type"
    )
    big_df = big_df.reindex(["pre", "post"], level="vaccination")

    grouped = DataFrame(
        data.groupby(["final_serotype_designation", "vaccination", "disease_type"])[
            "id"
        ].count()
    )

    big_df["count"] = grouped.id
    big_df.fillna(0, inplace=True)
    big_df.reset_index(inplace=True)

    # Adding a column to categorise serotypes into "plot" or "other"
    vaccine_serotypes = ["1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F"]
    plot_serotypes = list(data.final_serotype_designation.value_counts().head(10).index)

    for serotype in vaccine_serotypes:
        if serotype in plot_serotypes:
            pass
        else:
            plot_serotypes.append("{}".format(serotype))

    plot_status = []
    for x in big_df.index:
        if big_df.final_serotype_designation[x] in plot_serotypes:
            plot_status.append("plot")
        else:
            plot_status.append("other")

    big_df["plot_status"] = plot_status

    # pooling the "other" category
    other_df = big_df[big_df.plot_status == "other"]

    otherdf = DataFrame(
        other_df.groupby(["vaccination", "disease_type"])["count"].sum()
    )
    otherdf.reset_index(inplace=True)
    if data.country.unique()[0] == "Iceland":
        otherdf["sercat"] = [
            "other_pre",
            "other_pre",
            "other_pre",
            "other_pre",
            "other_post",
            "other_post",
            "other_post",
            "other_post",
        ]
    elif data.country.unique()[0] == "Kenya":
        otherdf["sercat"] = ["other_pre", "other_pre", "other_post", "other_post"]

    # joining the other category to the plotting df
    plot_df = big_df[big_df.plot_status == "plot"]
    plot_df.drop("plot_status", axis=1, inplace=True)

    plotdf = pd.concat([plot_df, otherdf])
    plotdf.reset_index(inplace=True)
    return plotdf


# Another function to generate the data required for stacked bars
def process_serostack(data):

    rawdata = serotype_process(data)

    subdflist = []
    for status in rawdata.sercat.unique():
        stack_count = []
        tempdf = rawdata[rawdata.sercat == status]
        rawvalues = list(tempdf["count"])

        if data.country.unique()[0] == "Iceland":
            # car
            stack_count.append(rawvalues[0])
            # inv
            stack_count.append(rawvalues[0] + rawvalues[1])
            # LRT
            stack_count.append(rawvalues[0] + rawvalues[1] + rawvalues[2])
            # OM
            stack_count.append(
                rawvalues[0] + rawvalues[1] + rawvalues[2] + rawvalues[3]
            )
        elif data.country.unique()[0] == "Kenya":
            stack_count.append(rawvalues[0])
            stack_count.append(rawvalues[0] + rawvalues[1])

        tempdf["stack_count"] = stack_count
        tempdf.sort_values(by=["stack_count"], ascending=False, inplace=True)

        subdflist.append(tempdf)

    # I think this huge concatenation is slowing everything down
    plot_data = pd.concat(subdflist, sort=True)
    plot_data.drop(
        ["count", "final_serotype_designation", "index", "vaccination"],
        axis=1,
        inplace=True,
    )
    return plot_data



##_____PLOTTING_____##

def stacked_serotype(data):

    # calling the functions and naming the plotting data here
    I_seroplot = process_serostack(data[data.country == "Iceland"])
    K_seroplot = process_serostack(data[data.country == "Kenya"])

    I_sub_seroplot = I_seroplot[I_seroplot.sercat != "other_pre"]
    I_sub_seroplot = I_sub_seroplot[I_sub_seroplot.sercat != "other_post"]

    K_sub_seroplot = K_seroplot[K_seroplot.sercat != "other_pre"]
    K_sub_seroplot = K_sub_seroplot[K_sub_seroplot.sercat != "other_post"]

    # function for sorting the data by overall frequency
    def serotype_order(data):

        vaccine_serotypes = ["1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F"]
        plot_serotypes = list(
            data.final_serotype_designation.value_counts().head(10).index
        )
        for serotype in vaccine_serotypes:
            if serotype in plot_serotypes:
                pass
            else:
                plot_serotypes.append("{}".format(serotype))

        serotype_order = []
        for serotype in data.final_serotype_designation.value_counts().index:
            if serotype in plot_serotypes:
                serotype_order.append("{}_pre".format(serotype))
                serotype_order.append("{}_post".format(serotype))
            else:
                pass

        return serotype_order

    # aesthetics
    sns.set_style("whitegrid")
    sns.set_context("paper")

    # Setting up axes and dimensions (2 rows, 1 col for landscape setup)
    stackedser_plot, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    plt.subplots_adjust(top=0.88, hspace=0.65, wspace=0.5)

    ax1data = I_sub_seroplot
    ax2data = K_sub_seroplot

    # setting up the order of the serotypes to be plotted
    I_order = serotype_order(data[data.country == "Iceland"])
    I_order.reverse()
    K_order = serotype_order(data[data.country == "Kenya"])
    K_order.reverse()

    # Plotting stacked bars
    sns.barplot(
        y="stack_count",
        x="sercat",
        hue="disease_type",
        data=ax1data,
        dodge=False,
        order=I_order,
        palette=[disease_pal[3], disease_pal[2], disease_pal[1], disease_pal[0]],
        ax=ax1,
    )
    sns.barplot(
        y="stack_count",
        x="sercat",
        hue="disease_type",
        data=ax2data,
        dodge=False,
        order=K_order,
        palette=[disease_pal[1], disease_pal[0]],
        ax=ax2,
    )

    # Title and legend tweaking
    ax1.set_title("Iceland", fontsize=14)
    ax2.set_title("Kenya", fontsize=14)
    for x in stackedser_plot.axes:
        x.set_ylabel("Pneumococci (n)", fontsize=12)
        x.set_yticklabels([0, 50, 100, 150, 200], fontsize=11, position=(0.01, 0))

    ax1.set_xlabel("")
    ax2.set_xlabel("Serotype (Pre/Post-PCV10 Introduction)", fontsize=12)

    ax1.set(xlim=[30, -1])
    ax2.set(xlim=[32, -1])

    ax1.get_legend().set_visible(False)
    ax2.get_legend().set_visible(False)

    # Manually defining a legend
    carbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[0], edgecolor="none")
    invbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[1], edgecolor="none")
    lrtbar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[2], edgecolor="none")
    ombar = plt.Rectangle((0, 0), 1, 1, fc=disease_pal[3], edgecolor="none")
    leg1 = ax1.legend(
        [carbar, invbar, lrtbar, ombar],
        ["Carriage", "Invasive Disease", "LRT Infection", "Otitis Media"],
        loc="upper right",
        ncol=4,
        prop={"size": 11},
        bbox_to_anchor=(1.02, -2.3),
    )
    leg1.draw_frame(True)

    # Adding vertical lines in between serotypes
    for x in stackedser_plot.axes:
        vlines = [1.5]
        linecount = 0
        if x is ax1:
            linecount = int(((len(ax1data.sercat.unique())) / 2) - 1)
        elif x is ax2:
            linecount = int(((len(ax2data.sercat.unique())) / 2) - 1)
        linecount_list = list(range(linecount))

        for value in linecount_list:
            vlines.append((2 * value) + 1.5)
        for value in vlines:
            x.axvline(value, 0, 1, color="black", linestyle="--", linewidth=1.2)

    # Adding serotype x axis labels
    for x in stackedser_plot.axes:
        x.xaxis.set_major_locator(ticker.IndexLocator(base=2, offset=0.8))
        serotypes = []

        if x is ax1:
            for sercat in I_order:
                ser = sercat.split("_")[0]
                # print(ser)
                if ser in serotypes:
                    pass
                else:
                    serotypes.append(ser)
        elif x is ax2:
            for sercat in K_order:
                ser = sercat.split("_")[0]
                # print(ser)
                if ser in serotypes:
                    pass
                else:
                    serotypes.append(ser)

        vaccine_serotypes = ["1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F"]
        style = []
        for ser in serotypes:
            if ser in vaccine_serotypes:
                style.append("bold")
            elif ser == "6E(6B)":
                style.append("bold")
            else:
                style.append("normal")

        x.set_xticklabels(
            serotypes, fontsize=11, position=(0, -0.15), va="top", rotation=90
        )

        for ser, st in zip(x.xaxis.get_ticklabels(), style):
            ser.set_fontweight(st)
        # x.tick_params(axis="y", direction="out", length=30, width=1, color="black")

    # Adding pre/post y axis labels
    for x in stackedser_plot.axes:
        data = []
        if x is ax1:
            data = ax1data
        elif x is ax2:
            data = ax2data
        length = len(data.sercat.unique())
        sublength = int(length / 2)
        vac_labels = sublength * ["post", "pre"]
        positionlist = list(range(length))
        positions = [0.6]
        for entry in positionlist:
            if entry == 0:
                pass
            else:
                positions.append(0.6 + entry)
        for vac, position in zip(vac_labels, positions):
            x.text(
                position,
                -5,
                "{}".format(vac),
                fontsize=10,
                rotation=70,
                verticalalignment="top",
            )

    # Adding zoomed in subplots
    ax3 = plt.axes([0.588, 0.68, 0.3, 0.18])  # this is [xcoord, ycoord, xsize, ysize]
    ax4 = plt.axes([0.653, 0.22, 0.235, 0.18])

    sns.barplot(
        y="stack_count",
        x="sercat",
        hue="disease_type",
        data=ax1data,
        dodge=False,
        order=I_order,
        palette=[disease_pal[3], disease_pal[2], disease_pal[1], disease_pal[0]],
        ax=ax3,
    )
    ax3.set(
        ylim=[0, 45],
        xlim=[11.5, -0.5],
        ylabel="",
        xlabel="",
        yticks=[0, 10, 20, 30, 40],
        xticklabels=[],
    )
    ax3.set_yticklabels(["", 10, 20, 30, 40], fontsize=12, position=(0.05, 0))

    sns.barplot(
        y="stack_count",
        x="sercat",
        hue="disease_type",
        data=ax2data,
        dodge=False,
        order=K_order,
        palette=[disease_pal[1], disease_pal[0]],
        ax=ax4,
    )
    ax4.set(
        ylim=[0, 55],
        xlim=[9.5, -0.5],
        ylabel="",
        xlabel="",
        yticks=[0, 10, 20, 30, 40, 50],
        xticklabels=[],
    )
    ax4.set_yticklabels(["", 10, 20, 30, 40, 50], fontsize=11, position=(0.05, 0))

    zoomplots = [ax3, ax4]
    for axis in zoomplots:
        axis.get_legend().set_visible(False)

    # Boxes around the plots
    for plot in stackedser_plot.axes:
        for _, spine in plot.spines.items():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1)

            if plot is ax3:
                spine.set_color("black")
                spine.set_linewidth(1)
            elif plot is ax4:
                spine.set_color("black")
                spine.set_linewidth(1)

    # Subplot horizontal lines between serotypes
    ax3_lines = [1.5, 3.5, 5.5, 7.55, 9.5]
    for x in ax3_lines:
        ax3.axvline(x, 0, 1, color="black", linestyle="--", linewidth=1.2)
    ax1.axvline(11.5, 0, 1, color="black", linestyle="-", linewidth=1)
    ax1.axvline(-0.5, 0, 1, color="black", linestyle="-", linewidth=1)

    ax4_lines = [1.5, 3.5, 5.5, 7.5]
    for x in ax4_lines:
        ax4.axvline(x, 0, 1, color="black", linestyle="--", linewidth=1.2)
    ax2.axvline(9.5, 0, 1, color="black", linestyle="-", linewidth=1)
    ax2.axvline(-0.5, 0, 1, color="black", linestyle="-", linewidth=1)

    return stackedser_plot


# function to describe how much of the overall dataset is not plotted in the stacked serotype figure 
# to be included in the figure legend whenever the figure is presented 
def seroplot_breakdown(data):
    def breakdown(data):

        vaccine_serotypes = ["1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F"]
        plot_serotypes = list(
            data.final_serotype_designation.value_counts().head(10).index
        )
        for serotype in vaccine_serotypes:
            if serotype in plot_serotypes:
                pass
            elif serotype in data.final_serotype_designation.value_counts().index:
                plot_serotypes.append("{}".format(serotype))
            else:
                pass

        data.final_serotype_designation.value_counts().head(10)
        plot_number = sum(data.final_serotype_designation.value_counts().head(10))

        print("Number of serotypes included in plot:")
        print(len(plot_serotypes))

        print("Number of isolates represented in plot:")
        print(plot_number)

        print("Number of serotypes not included in plot:")
        print(len(data.final_serotype_designation.value_counts()) - len(plot_serotypes))

        print("Number of isolates not represented in plot:")
        print(len(data) - plot_number)

    print("Iceland:")
    breakdown(data[data.country == "Iceland"])
    print("____\n")
    print("Kenya:")
    breakdown(data[data.country == "Kenya"])



