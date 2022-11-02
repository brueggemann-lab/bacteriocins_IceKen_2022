###
#####_____CONTIGUITY CAT_____#####
###

# bin each cluster based on contiguity
# for each cluster present (not fragment) pull out the expected loci
# for each locus, pull out the contig and the start and stop location
# categorise according to the following bins

# contiguous
# all loci are on the same contigs
# within Xbp of the expected adjacent locus
# tolerated separation of loci is determined by threshold argument, deafult = 2500bp 

# EOC 
# the cluster is found on multiple contigs 
# is in multiple fragments which are themselves contiguous 
# fragments are found within Xbp of the contig break 
# tolerated distance from EOC determined by the eoc_threshold argument, default = 1500bp

# contiguous with Ns 
# the loci coordinates suggest contiguity
# but some ambiguous reads so the coordinates do not necessarily reflect "real" contiguity 
# these are counted as contiguous 

# non-contiguous 
# various sub-categories associated with this one 
# depending on same contig vs. multiple contigs
# and whether or not they are EOC-adjacent 

# final output: summary table describing the frequency of each contiguity category for each bacteriocin
# the summary is useful for supp. tables describing the result
# also a full table showing the results for each isolate 
# the full export is used for when the cluster contiguity is incorporated in analysis code 

#####____________set up - library imports____________#####

from click.decorators import pass_context
import pandas as pd
from pandas import DataFrame, Series
import numpy as np

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation


#####____________set up - command line interface____________#####
# click
# bigs_export is a dataset export of all the isolates which the contiguity check is checking
# including curated bacteriocin loci

@click.command()
# file name of BIGSdb export 
@click.option(
    "--bigs_export",
    required=True,
    help="file exported from BIGSdb with curated bacteriocin loci",
)
# path to the directory which contains the gbk-like files with locus coordinates
# important that this has a "/" at the end  
@click.option(
    "--gbk_dir",
    required=True,
    help="location of the annotated contig files in gbk format",
)

@click.option(
    "--threshold",
    default=2500,
    help="upper limit of proximity of loci to count as contiguous",
)
@click.option(
    "--eoc_threshold",
    default=2500,
    help="maximum distance from EOC for a locus for that cluster to be considered ambiguous",
)

###
#####____________main()____________#####
###

def main(bigs_export, gbk_dir, threshold, eoc_threshold):
    #  calls the contiguity_cat function
    # imports the bigs export
    # processes it
    # performs contiguity checks
    # returns a summary for contguity
    contiguity_cat(bigs_export, gbk_dir, threshold, eoc_threshold)
    return

###
#####____________master function calling all the others for the given data export____________#####
###

def contiguity_cat(bigs_export, gbk_dir, threshold, eoc_threshold):

    # process the dataframe
    # using functions defined below
    df = process(bigs_export)

    # categorise the clusters based on contiguity 
    full_cat_df = categorise(df, gbk_dir, threshold, eoc_threshold)
    # for the full output, only taking the cluster contiguity cat columns 
    # to be used in Narwhal 
    cat_cols = [col for col in full_cat_df if col.endswith("contiguitycat")]
    cat_cols.append("id")
    slice_cat_df = full_cat_df[cat_cols]

    # save the full output
    output_full_name = "cluster_cont_" + bigs_export.split("_")[1] + "_" + str(threshold) + "_" + str(eoc_threshold) + ".csv"
    output_full_path = "./outputs/" + output_full_name 
    slice_cat_df.to_csv(output_full_path)

    # summarise the categories 
    output_summary = contiguity_summary(full_cat_df)

    # logging
    print("Generated summary of contiguity in bacteriocin clusters:")
    
    # print the summary in command line interface 
    print(output_summary)

    # save the summary of cluster contiguity 
    # for use as supplementary data 
    output_summary_name = "cluster_cont_" + bigs_export.split("_")[1] + "_" + str(threshold) + "_" + str(eoc_threshold) + "_summary.csv"
    output_summary_path = "./outputs/" + output_summary_name 
    output_summary.to_csv(output_summary_path)

###
#####____________importing the bigs data export____________#####
###

def data_import(bigs_export):
    input_file = "./" + bigs_export
    # logging
    print("Reading in file: " + input_file)
    return pd.read_csv(input_file, low_memory=False)


###
#####____________processing the dataframe____________#####
###

# processing taken from Narwhal

# making sure the allele designations are strings
# essential for the definition of the bacteriocin profiles
def df_types(data):
    for x in data.columns:
        if x.startswith("mj"):
            data["{}".format(x)] = data["{}".format(x)].astype(str)
        else:
            pass
    # logging
    print("Adjusted allele designations to string type")

# specifying the order of the loci in the clusters
def adjust_order(data):

    # first specify the order of the loci
    locus_order = []

    # dict with the order of the non-alphabetic loci
    # does not include bacteriocins where loci order is alphabetical (e.g. streptococcins, ABC)

    non_alph = {
        "sla": [
            "mj_slaA1",
            "mj_slaA2",
            "mj_slaA3",
            "mj_slaA4",
            "mj_slaA5",
            "mj_slaF",
            "mj_slaE",
            "mj_slaK",
            "mj_slaR",
            "mj_slaM",
            "mj_slaT",
        ],
        "slb": ["mj_slbF", "mj_slbG", "mj_slbE", "mj_slbA", "mj_slbM", "mj_slbT"],
        "slc": ["mj_slcA", "mj_slcX", "mj_slcL", "mj_slcT"],
        "sle": [
            "mj_sleM1",
            "mj_sleA1",
            "mj_sleA2",
            "mj_sleM2",
            "mj_sleM3",
            "mj_sleT",
            "mj_sleX1",
            "mj_sleF",
            "mj_sleG",
            "mj_sleX2",
        ],
        "slg": [
            "mj_slgA1",
            "mj_slgA2",
            "mj_slgM",
            "mj_slgD",
            "mj_slgP1",
            "mj_slgT",
            "mj_slgP2",
        ],
        "slh": [
            "mj_slhP",
            "mj_slhR",
            "mj_slhK",
            "mj_slhF",
            "mj_slhE",
            "mj_slhG",
            "mj_slhX1",
            "mj_slhX2",
            "mj_slhA",
            "mj_slhB",
            "mj_slhT",
            "mj_slhC",
            "mj_slhI",
        ],
        "sli": [
            "mj_sliP",
            "mj_sliR",
            "mj_sliK",
            "mj_sliF",
            "mj_sliE",
            "mj_sliG",
            "mj_sliA",
            "mj_sliB",
            "mj_sliT",
            "mj_sliC",
            "mj_sliI",
        ],
        "slj": [
            "mj_sljA1",
            "mj_sljL",
            "mj_sljP",
            "mj_sljT1",
            "mj_sljT2",
            "mj_sljT3",
            "mj_sljA2",
        ],
        "sls": [
            "mj_slsA",
            "mj_slsC",
            "mj_slsB1",
            "mj_slsB2",
            "mj_slsF",
            "mj_slsE",
            "mj_slsG",
            "mj_slsR",
            "mj_slsK",
        ],
        "ssa": ["mj_ssaA", "mj_ssaCD", "mj_ssaX1", "mj_ssaX2", "mj_ssaP", "mj_ssaX3"],
    }

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


# changing the order of the columns on the input dataframe
# called in process()
def col_order(data):
    data = data[adjust_order(data)]
    print(
        "Adjusted locus order in curated data export to reflect order found in refernce clusters"
    )
    return data

# functions returning the unique clusters and loci in the dataframe
# may not be needed, remove if not used
# taken from Narwhal, where they are useful
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
    # this function can be modified if the loci naming convention in BIGS is changed
    loci = [x for x in data.columns if x.startswith("mj")]
    return loci

# bacteriocin profile class
# taken from Narwhal, simplified for use here
# updated to meet class naming conventions 
class BacteriocinProfile:
    def __init__(self, dict):

        self.cluster = list(dict.keys())[0][:3]

        self.allelic_profile = "-".join(dict.values())

        self.partfull_profile = "-".join(
            ["/" if v == "0" else k[3:] for k, v in dict.items()]
        )

        self.category = cluster_category(self.partfull_profile)

        self.status = True if self.category in ["P", "Fl"] else False


# allele_dict generates a dictionary with loci and alleles for a given cluster from the dataframe of curated data
# for input into the BacteriocinProfile class
def allele_dict(data, cluster, index):

    loci = [x[3:] for x in data.columns if cluster in x]
    output = {}
    for col in loci:
        output[col] = data["mj_{}".format(col)][index]

    return output

# function for generating the generalised full/partial status attribute in the BacteriocinProfile class
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

# this returns a dataframe with the index of the dataframe with columns for all the BacteriocinProfile attributes
# this can be merged onto the dataframe
def cluster_cols(data, cluster):

    cluster_df = DataFrame(
        index=data.index,
        columns=[
            "{}_a_profile".format(cluster),
            "{}_profile".format(cluster),
            "{}_category".format(cluster),
            "{}_status".format(cluster),
        ],
    )

    for x in data.index:
        # print(x)
        cluster_df["{}_a_profile".format(cluster)][x] = BacteriocinProfile(
            allele_dict(data, cluster, x)
        ).allelic_profile
        cluster_df["{}_profile".format(cluster)][x] = BacteriocinProfile(
            allele_dict(data, cluster, x)
        ).partfull_profile
        cluster_df["{}_category".format(cluster)][x] = BacteriocinProfile(
            allele_dict(data, cluster, x)
        ).category
        cluster_df["{}_status".format(cluster)][x] = BacteriocinProfile(
            allele_dict(data, cluster, x)
        ).status

    return cluster_df

# this function calls cluster_cols for all the bacteriocins of interest
# and merges the new cols onto the given dataframe
# returns a separate object
def profile_cols(data):

    output = data
    for x in unique_clusters(data):
        print(x)
        output = output.merge(
            cluster_cols(data, x), how="outer", left_index=True, right_index=True
        )
    # logging
    print("Added cluster columns")

    return output

# processing function taken directly from Narwhal.ipynb
# deleted parts that are not required here
# added the bigs_export arg
# inherited from the command line arg
# path to the bigs_export csv
def process(bigs_export):
    # imports the dataframe
    rawdf = data_import(bigs_export)

    # adjusting types
    df_types(rawdf)
    # adjusting order of the locus columns
    rawdf = col_order(rawdf)
    # adding columns with bacteriocin profiles using the profile class attributes
    rawdf = profile_cols(rawdf)

    print("Processed bigs data export, ready for contiguity checks")
    return rawdf


###
#####____________CONTIGUITY CHECKS____________#####
###

#####____________useful functions____________#####

# function to return the expected loci of a cluster based on the processed dataframe
# called based on index
# loci returned in the cluster order, not alphabetical
# assuming forward direction
def expected_loci(index, cluster, data):
    expected_loci = []
    for x in data["{}_profile".format(cluster)][index].split("-"):
        if x == "/":
            pass
        else:
            expected_loci.append(("mj_{}".format(cluster) + x))
    # [ for x in df["{}_profile".format(cluster)][index].split("-")]
    return expected_loci


#####____________fetching locus coordinates____________#####

# function to get all the coordinates for the loci of a given cluster
# returns a dict describing the coordinates of the loci of a single cluster
# in a single isolate
# only on the loci that are detected in that genome according to the curation

def coordinates(index, cluster, seq_records, data):
    # empty dict, keys will be the expected loci, populated below
    # each entry is another dict with the contig number and the start and end point of the locus
    coords = {}
    # is being fed each cluster
    # expected loci should return the correct loci
    # find the annotations in the genbank file
    for record in seq_records:
        for feature in record.features:
            locus = feature.qualifiers["locus"][0]

            if locus in expected_loci(index, cluster, data):

                # writing an exception for where a locus could not be parsed
                # these are cases where the coordinates extend before the start of the contig
                # i.e. negative coordinates
                # in cases where the locus is 5' EOC

                if str(feature.location) == "None":
                    # give the locus coordinates of 0
                    coords[locus] = {
                        "contig": record.id,
                        "start": 0,
                        "end": 0,
                        "direction": 0,
                    }

                else:
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    direction = feature.strand
                    coords[locus] = {
                        "contig": record.id,
                        "start": start,
                        "end": end,
                        "direction": direction,
                    }
    return coords


#####____________defining contiguity check functions____________#####

# a set of functions which use the coordinates of the loci to determine the contiguity category
# returns the category
# thresholds set here as global variables, could be set as arguments with CLICK

# threshold = 2500
# eoc_threshold = 2500

# generates a dict desribing which contig each locus is found on 
# based on the coords output of the coordinates function 
def contig_dict(coords):
    output = {}
    for locus in coords.keys():
        if coords[locus]["contig"] in output.keys():
            if locus in output[coords[locus]["contig"]]:
                pass
            else:
                output[coords[locus]["contig"]].append(locus)
        else: 
            output[coords[locus]["contig"]] = [locus]
    return output

# contiguity check works on a single cluster in a single isolate
# isolate accessed through the index
def contiguity_check(index, cluster, seq_records, data, threshold, eoc_threshold):
    
    # first of all we don't care about a cluster that is absent
    # nor fragment clusters 
    # so use the cluster_status column to find these 
    if data["{}_status".format(cluster)][index]:
        pass
    else:
        return "Absent"
    
    # get the seq coordinates and contig number
    # of the expected loci from the given cluster
    # in the given isolate
    coords = coordinates(
        index, cluster, seq_records, data
    )  # **feed this the parsed seqrecords

    # first checking that the loci are on the same contig
    # each contig with a list of the loci of the cluster found on that contig 
    contigs = contig_dict(coords)

    if len(contigs.keys()) == 1:
        # check for adjacency
        # returns whatever adjacency returns
        return adjacency(index, cluster, coords, threshold, data)
    else:
        # check for proximity to the end of the contig
        return end_of_contig(coords, seq_records, threshold, eoc_threshold)
    # need to adjust the end of contig function to check for all contigs 

# adjacency uses a set of coords for the loci of a cluster
# to determine whether the loci are found adjacent to one another
# assumes loci are on the same contig
# returns categories accordingly
def adjacency(index, cluster, coords, threshold, data):
    loci = expected_loci(index, cluster, data)
    
    # then go on to check for adjacency 
    for position, locus in enumerate(loci):
        # skip the last locus
        if position == (len(loci)) - 1:
            pass
        # checking the start and end positions of each locus
        # in the order they are found in a cluster
        else:
            # pull out end of first locus
            locus1_end = coords[loci[position]]["end"]
            # pull out start of adjacent locus
            locus2_start = coords[loci[position + 1]]["start"]

            # exception handling reverse orientation
            if coords[locus]["direction"] == -1:
                locus1_end = coords[loci[position]]["start"]
                locus2_start = coords[loci[position + 1]]["end"]

            # find the difference
            # abs() forces positive
            # takes into account the possibility of reverse orientation
            difference = abs(locus2_start - locus1_end)
            # checking proximity
            if difference > threshold:
                return "Non-contiguous (one contig)"
            else:
                pass
    
    # if any locus has the "N" flag
    # needs to be binned to an ambigous read category       
    # scaffolding with Ns may prevent certainty in contiguity 
    for locus in loci: 
        if data["{}".format(locus)][index] == "[S]":
            return "Contiguous with Ns"
        else:
            pass
    # if none of the pairs of loci are non-contiguous,
    # and no ambiguous reads 
    # the cluster is counted as contiguous
    return "Contiguous"

def end_of_contig(coords, seq_records, threshold, eoc_threshold):      #### TO DO: ADD CHECK THAT LOCI ARE CONTIGUOUS EITHER SIDE OF CONTIG BREAK #####
    # seq_records arg is the parsed gbk file for this isolate

    # need to check whether any locus is within X kbp of EOC
    # if so, cluster may be contiguous but we can't tell, call it "EOC" 
    # if in the middle of multiple contigs, non-contiguous
    # i.e. all loci > X kbp away from EOC
    contigs = contig_dict(coords)
    
    # list of contigs which have an EOC locus 
    # should be the same as contigs.keys() if genuine EOC category 
    EOC_contigs = []

    # look at each contig one at a time 
    for contig in contigs.keys():
        # placeholder for contig length
        contig_length = 0 
        for record in seq_records:
            if record.id == contig:
                contig_length = len(record.seq)
        
        # check that contig length is found
        if contig_length == 0:
            print("error: can't find contig length")
            # if length is still 0, something has gone wrong
            return "error"
        else:
            pass

        # add to this each time a locus close to EOC is observed
        # close meaning within the EOC threshold 
        EOC_count = 0
        # list of loci on the contig 
        loci = contigs[contig]
        # look at position of each locus in that contig 
        for locus in loci: 
            if coords[locus]["start"] < eoc_threshold:
                EOC_count = EOC_count + 1 
            elif coords[locus]["start"] > (contig_length - eoc_threshold):
                EOC_count = EOC_count + 1
            elif coords[locus]["end"] < eoc_threshold:
                EOC_count = EOC_count + 1
            elif coords[locus]["end"] > (contig_length - eoc_threshold):
                EOC_count = EOC_count + 1
            else:
                pass
        if EOC_count >= 1:
            EOC_contigs.append(contig)
        else:
            pass
    
    ## need to change how this works 
    # check that all loci at each EOC are within the contiguity threshold 

    # use EOC_contigs - populated list of contigs which have a locus at the end 
    # go through the list
    # get all loci off a single contig and do a contiguity check - like adjacency()

    
    # coords - dict with keys for each locus, each entry also a dict 
    if len(EOC_contigs) == len(list(contigs.keys())):
        
        for contig in EOC_contigs: 

            # pull out loci on this contig 
            # check adjacency to each other 
            
            # contig_loci is a sub-section of coords
            # only containing coordinates of loci on the contig 
            contig_loci_coords = {}
            for locus in coords.keys():
                if coords[locus]["contig"] == contig:
                    contig_loci_coords[locus] = (coords[locus]) 
                else:
                    pass
            
            # all of this is adapted from adjacency() 
            # list of loci is contig_loci.keys() 
            contig_loci = list(contig_loci_coords.keys())
            for position, locus in enumerate(contig_loci):
                # skip the last locus
                if position == (len(contig_loci)) - 1:
                    pass
                else:
                # checking the start and end positions of each locus
                # in the order they are found in a cluster
                # pull out end of first locus
                    locus1_end = contig_loci_coords[contig_loci[position]]["end"]
                # pull out start of adjacent locus
                    locus2_start = contig_loci_coords[contig_loci[position + 1]]["start"]
                # exception handling reverse orientation
                    if contig_loci_coords[locus]["direction"] == -1:
                        locus1_end = contig_loci_coords[contig_loci[position]]["start"]
                        locus2_start = contig_loci_coords[contig_loci[position + 1]]["end"]
                # find the difference
                # takes into account the possibility of reverse orientation
                    difference = abs(locus2_start - locus1_end)
                    # checking proximity
                    if difference > threshold:
                        return "Non-contiguous (multiple contigs, non-adjacent loci)"
                    else:
                        pass
        # return the "EOC" category after the adjacency of loci on each contig has been checked 
        return "EOC"
    
    else:
        # returns this if the fragments are on multiple contigs and not near the contig break 
        return "Non-contiguous (multiple contigs, not EOC-adjacent)"            

#####____________performing contiguity checks____________#####

### function calling above functions and returning a dataframe
# writes the outputs of the above functions into an empty dataframe
#  merges the category dataframe onto the processed df export


def categorise(data, gbk_path, threshold, eoc_threshold):
    cat_cols = [(x + "_contiguitycat") for x in unique_clusters(data)]

    output_df = DataFrame(index=data.index, columns=cat_cols)

    for index in data.index:
        # parse the relevant gbk file once
        path = gbk_path + "id_" + str(data.id[index]) + "_contigs.gbk"
        seq_records_raw = SeqIO.parse(path, "genbank")
        seq_records = list(seq_records_raw)
        print("Parsed seq record for isolate " + str(data.id[index]))
        # access the parsed gbk file multiple times
        for cluster in unique_clusters(data):

            ## this is passing contiguity_check the right cluster
            output_df["{}_contiguitycat".format(cluster)][index] = contiguity_check(
                index, cluster, seq_records, data, threshold, eoc_threshold
            )

    # returns the processed dataframe with the category columns added on
    return pd.concat([data, output_df], axis=1)


#####____________summarising cluster contiguity____________#####

#  generating a summary of the contiguity categories
# takes the processed dataframe with the category columns as the input
def contiguity_summary(cat_data):
    # populates an empty list with dfs
    # each one is a summary of a cluster
    df_list = []
    for cluster in unique_clusters(cat_data):
        df_list.append(cluster_summary(cat_data, cluster))
    # merge each cluster df into a single summary df
    output = pd.concat(df_list)
    return output


# generates a summary for a single cluster
# to be run on a merged dataframe with the cat columns
def cluster_summary(catdf, cluster):
    vc = DataFrame(catdf["{}_contiguitycat".format(cluster)].value_counts())

    if "Absent" in vc.index:
        if len(vc.index) == 1:
            # cluster was not detected, nothing to count
            return
        else:
            # remove absent, not interested
            vc = vc[vc.index != "Absent"]
    else:
        # do nothing if cluster is never absent
        pass

    # calculating the % of observed clusters that fall into each category
    total = sum(vc["{}_contiguitycat".format(cluster)])
    perc_col = Series(index=vc.index, dtype=np.float64)
    for status in vc.index:
        count = vc["{}_contiguitycat".format(cluster)][status]
        perc_col[status] = 100 * ((count) / total)

    # rename columns
    vc.columns = ["Frequency"]
    vc["Proportion_of_clusters"] = perc_col

    # generating a multiindex with both cluster name and contiguity category
    index = []
    for x in vc.index:
        index.append(("{}".format(cluster), x))
    multiindex = pd.MultiIndex.from_tuples(index, names=["Cluster", "Category"])

    vc = vc.set_index(multiindex)
    return vc


#####____________calling main()____________#####

main()
