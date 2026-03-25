## IMPORT MODULES ##

import numpy as np
import os  # for saving to folder
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcdefaults() # default settings to avoid warnings related to font.

import socket
import os as os
import sys as sys
import multiprocessing as mp
from pysam import AlignmentFile
import argparse

import requests # to look for weblinks. 
from bs4 import BeautifulSoup # to transform into table information in a database.
import urllib # also to retrieve weblinks.
import time
from liftover import get_lifter # to convert positions from hg38 to hg19 (used to align reads).
import seaborn as sns # for visualization.
import re # to work with regular expressions.
import pathlib
from pathlib import Path

## IMPORT FUNCTIONS ##

from pulldown import load_snp_file_OY, ref_alt_count, pulldown_bamtable, call_y_bam, create_tree, div_anc_der, create_parent_dct, create_path, network, unique_lineages

def main():

    # Define parameters for the algorithm.  
    parser = argparse.ArgumentParser(
        description = "Assign most derived Y haplogroups in samples."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help = "Provide a TSV file with valid routes to .bam files [columns: iid, bam, sex (f or m)]."
    )

    parser.add_argument(
        "--index_i",
        required=False,
        help = "Specify which is the first sample in the set to analyse [default 0]."
    )

    parser.add_argument(
        "--index_f",
        required=False,
        help = "Specify which is the last sample in the set to analyse [defult len(dataset)]."
    )

    parser.add_argument(
        "--base_qual",
        required=False,
        help = "Minimum base quality [default 20]"
    )

    parser.add_argument(
        "--map_qual",
        required=False,
        help = "Minimum mapping quality [default 25]"
    )

    parser.add_argument(
       "--snps",
       required=False,
       help = "Provide a .csv file containing information for SNPs in Y chromosome [pos,ref,alt,ID]"
    )

    parser.add_argument(
        "--rg",
        required=False,
        help = "Specify the reference genome for genome coordinates (either hg38 or hg19) [default = hg38]."
    )

    parser.add_argument(
        "--create_network",
        required=True,
        help = "Create a network connecting the different Y-haplogroups and placingsamples in it [Y/N]"
    )

    parser.add_argument(
        "--width",
        required=False,
        help = "Width of the network created [default 120]."
    )

    parser.add_argument(
        "--height",
        required=False,
        help = "Height of the network created [default 90]."
    )


    # Store arguments as variables:
    args = parser.parse_args()

    bam = args.bam

    # If no set of SNPs provided, use the default dataset.
    path_snps = args.snps
    if path_snps is None:
        path_snps = "data/input/YBrowse_snps_hg38-runtime.csv"

    # Only when using default set of markers, consider genomic positions on hg38 reference.
    reference_genome = args.rg
    if reference_genome is None:
        reference_genome = "hg38"

    if args.index_i is not None: # only if a starting index is provided.
        initial = args.index_i
    else:
        initial = 0 # if not, start the analysis in the first sample.

    if args.index_f is not None:
        final = args.index_f
    else:
        final = 0 # if no final is provided, store 0 for the object to then sum the total number of samples (i.e. analyse everything).

    base_qual = args.base_qual
    map_qual = args.map_qual
    create_network = args.create_network
    width = args.width
    height = args.height

    print("\n--- Running y_call software ---\n ")

    # Create a dictionary to store relations between child (key) and parent (value) haplogroups.
    chpar = create_parent_dct()

    # Path to the set of markers considered in the analysis (after filtering). 
    file_OY = Path("./data/input/all_snps.csv")

    if file_OY.exists():
        print("Dataset for filtered markers already saved: ./data/input/all_snps.csv")
        OY = pd.read_csv(file_OY)
        path_bed = Path("./data/output/OY_snps.bed")
        if path_bed.exists:
            print(f"Corresponding .bed file in {path_bed}")
        else:
            savepath = "./data/output/OY_snps.bed" # save a BED file with coordinates for every SNP. 

            dft = OY[["chrom", "pos"]].copy()
            dft["pos1"] = dft["pos"]
            dft.to_csv(savepath, sep="\t", index=False, header=None)
            print(f"Coresponding .bed file in {path_bed}")
    else:
        print("Dataset for filtered markers missing\nCreating file...")
        OY = load_snp_file_OY(path_snps, reference_genome, chpar, unique=True) # import all SNPs in .csv file. 
        OY.to_csv("data/input/all_snps.csv")
        print(f"Done\nSaved as ./data/input/all_snps.csv")

        OY = OY.sort_values(by="pos") # sort by position.

        savepath = "./data/output/OY_snps.bed" # save a BED file with coordinates for every SNP. 

        dft = OY[["chrom", "pos"]].copy()
        dft["pos1"] = dft["pos"]
        dft.to_csv(savepath, sep="\t", index=False, header=None)
        print(f"Coresponding .bed file in {savepath}")

    dft = pd.read_csv(bam, sep="\t") # file with all routes to BAM files.
    
    n = dft[dft["bam"].notna()] # select only rows where BAM file is available. 
    males = n[n["sex"] == "m"] # select only male samples.
    
    bam_dict_males = dict(zip(males["iid"], males["bam"])) # create a dictionary as [SAMPLE: route .bam].
    bam_dict_males["PTN261_ss.A0101"] = "/mnt/archgen/Autorun_eager/eager_outputs/SG/PTN/PTN261/merged_bams/initial/PTN261_ss_udghalf_libmerged_rmdup.bam"
    
    if final == 0:    
        routes = list(bam_dict_males.items())[int(initial):int(final+len(bam_dict_males)+1)] # if no indexes provided, consider all samples. 
    else:
        routes = list(bam_dict_males.items())[int(initial):int(final)] # if indexes provided, select only those samples.     

    print(f"\n- Loaded a total of {len(routes)} male individuals. -\n")

    # Store in a list content for the final dataFrame:
    samples = []
    branches = []
    levels = []
    ratios = []
    scores = []

    # Create a dictionary to store the total number of SNPs per haplogroup:
    dict_snps = {}

    for sample in routes:

        # Show content for the sample being analysed.
        iid = sample[0]
        bam = sample[1]

        if not Path(bam).exists(): # only if .bam file was not found.
            print(f"Sample {sample} not consdiered, missing BAM file.")
            continue

        print(f"\n- Summary coverage statistics for sample {iid} -")

        # Do the call, depending on if the user has specified a minimum base quality and mapping quality.
        if base_qual is not None and map_qual is not None:
            df_ch, df_der = call_y_bam(df=OY, path_bam=bam,
                                       path_bed='data/output/OY_snps.bed',
                                       path_temp=f'data/output/temp/{iid}_temp.tsv',
                                       snip5=0, snip3=0, base_qual=base_qual, map_qual=map_qual)

        elif base_qual is not None:
            df_ch, df_der = call_y_bam(df=OY, path_bam=bam,
                                       path_bed='data/output/OY_snps.bed',
                                       path_temp=f'data/output/temp/{iid}_temp.tsv',
                                       snip5=0, snip3=0, base_qual=base_qual)

        elif map_qual is not None:
            df_ch, df_der = call_y_bam(df=OY, path_bam=bam,
                                       path_bed='data/output/OY_snps.bed',
                                       path_temp=f'data/output/temp/{iid}_temp.tsv',
                                       snip5=0, snip3=0, map_qual=map_qual)
        else:
            df_ch, df_der = call_y_bam(df=OY, path_bam=bam,
                                       path_bed='data/output/OY_snps.bed',
                                       path_temp = f'data/output/temp/{iid}_temp.tsv',
                                       snip5=0, snip3=0)

        # Post-process filtering (subset of SNPs to be exluded in data/input/mm_v3.tsv)
        df_cts = pd.read_csv("data/input/mm_v3.tsv", sep="\t")
        df_ex3 = df_cts[df_cts["count"]>5] # SNPs that are ancestral in at least 5/209 test sample derived chains.

        # Create a dataFrame to know the total number of SNPs per branch, and the corresponding ancestral and derived counts. 
        dft = div_anc_der(df_ch, df_exclude=df_ex3)
        dfd = dft 

        # Add a new column to the dataset as the ratio of derived SNPs in respect to the Total SNPs.
        dfd["Derived/Total"] = round(dfd["Derived"] / dfd["Total_SNPs"],6)

        # Add a new column to the dataset as the ratio of #ANC in par. in #DER in par.
        dfd["Ratio"] = round(dfd["#ANC in par."] / dfd["#DER in par."],6)

        # Add a new column to the dataset as the score for every haplogroup = #DER in par. + Derived - 3* (#ANC in par. + Ancestral).
        dfd["Score"] = round((dfd["#ANC in par."]*-3)+(dfd["#DER in par."])+(dfd["Ancestral"]*-3)+(dfd["Derived"]),6)

        # Exclude haplogroup I-P38, as it's got a different name.  
        df = dfd[dfd["Branch"]!="I-P38"].sort_values(by="Score")
        #display(df.tail(5))

        # Create directories for SNP and Y-haplogrooups information.
        file_path = pathlib.Path(f"data/output/markers/{iid}")
        file_path.mkdir(parents=True, exist_ok=True)

        file_path2 = pathlib.Path(f"data/output/haplogroups/{iid}")
        file_path2.mkdir(parents=True, exist_ok=True)

        # Write a .tsv file to consult SNPs for every branch.
        df_ch.to_csv(f"data/output/markers/{iid}/snps_{iid}.tsv", sep="\t")
        print(f"\nSaved data for pileup as data/output/markers/{iid}/snps_{iid}.tsv")

        # Write a .tsv file to consult only derived SNPs for every branch.
        df_der.to_csv(f"data/output/markers/{iid}/snps_der_{iid}.tsv", sep="\t") 
        print(f"Saved data for pileup (considering only derived sites) as data/output/markers/{iid}/snps_der_{iid}.tsv") 

        # Write a .tsv file to consult information for every Y-haplogroup called. 
        df.to_csv(f"data/output/haplogroups/{iid}/haplogroup_{iid}.tsv", sep="\t")

        # Write entries in the dictionary to see the total of SNPs per haplogroup.
        for haplo, total_snps in zip(df["Branch"], df["Total_SNPs"]):
            if haplo in dict_snps:
                dict_snps[haplo].append(int(total_snps))
            else:
                dict_snps[haplo] = [int(total_snps)]

        # Create directory for the tree files created in each sample.
        file_path3 = pathlib.Path(f"data/output/trees/{iid}")
        file_path3.mkdir(parents=True, exist_ok=True)

        # Create a tree to know where is the sample located in Y-haplogroups
        with open(f"data/output/trees/{iid}/tree_output_{iid}.txt", "a") as f:
            f.write(f"Tree created for sample {iid}, with most derived Y-haplogroup: {df.tail(1)["Branch"].iloc[0]}\n")
            create_tree(string = df.tail(1)["Branch"].iloc[0],level = int(df.tail(1)["Level"].iloc[0]), df=df, chpar=chpar, file=f)
            f.write("\n")

        # Create a file to store all paths (for every sample), useful to then create the final tree
        with open(f"data/output/paths.txt", "a") as f2:
            path = []
            create_path(string = df.tail(1)["Branch"].iloc[0], df=df, chpar=chpar, path=path)
            f2.write(iid.split("_")[0] + ": " + ",".join(path))
            f2.write("\n")
        print(f"Y-haplogroup simplified tree saved as data/output/trees/{iid}/tree_output_{iid}.txt\nPath to most derived haplogroup included in data/output/paths.txt")

        # Keep info for summary
        samples.append(iid)
        branches.append(df.tail(1)["Branch"].iloc[0])
        levels.append(int(df.tail(1)["Level"].iloc[0]))
        ratios.append(float(df.tail(1)["Ratio"].iloc[0]))
        scores.append(int(df.tail(1)["Score"].iloc[0]))

    # Create a final dataFrame to store all measures for every sample.
    summary_df = pd.DataFrame({
        "Sample": samples,
        "Most derived": branches,
        "Level": levels,
        "Ratio": ratios,
        "Score": scores
    })

    # Write a .csv file with the summary.
    summary_df["Sample"] = summary_df["Sample"].str.split("_").str[0]
    summary_df.to_csv("data/output/scores.csv", mode="a", header=False, index=False)
    print(f"\nSaved Y-haplogroup calls for {len(samples)} individuals as data/output/scores.csv.")

    # Write a .csv file with SNPs per branch.
    data_list = []
    for haplo, snps in dict_snps.items():
        for snp in snps:
            data_list.append({"Haplogroup": haplo, "Total_SNPs": snp})

    # Create DataFrame and save.
    df_snps = pd.DataFrame(data_list)
    df_snps.to_csv("data/output/snps_branch.csv", mode="a", header=False, index=False)

    # Compute the total number of unique lineages for the dataset.
    # Use info for different paths in each sample and the corresponding scores. 
    print("\n- Running analysis to know the total number of unique lineages in dataset -\n")
    with open ("data/output/paths.txt", "r") as file:
        items = file.readlines()

    data = []
    for i in items:
        data.append(i.strip())

    summary_df = pd.read_csv("data/output/scores.csv",header=None)
    summary_df.columns = ["Sample","Most derived","Level","Prop.","Score"]

    unique = unique_lineages(summary_df, data)

    print(f"\nTotal number of unique lineages in dataset: {len(unique)}.\n")

    # Create a network to place samples in a Y-haplogroup tree (only if specified by the user).
    if create_network == "Y":
        print("\n- Creating a network to know relationship between samples. -\n")
        df = pd.read_csv("data/output/snps_branch.csv", header=None)

        # Divide the dataFrame by the different haplogroups.
        df.columns = ["Haplogroup", "SNPs_num"]
        grouped = df.groupby("Haplogroup")

        # Generate a new dataFarme with columns for the Branch, the Level and the average total of SNPs.
        new_df = grouped.agg(
                Sum=("SNPs_num", "mean")
        ).round(2).reset_index()

        avg_snps = dict(zip(new_df["Haplogroup"], new_df["Sum"]))

        if width is not None and height is not None:
            network(data, int(width), int(height), avg_snps)
        else:
            network(data, 120, 90, avg_snps) # if imensions not provided by the user, use default options.
        print(f"\nTree created for a total of {len(data)} samples\n")


if __name__ == "__main__":

    """ Run function main() only when the program is called from the CLI"""

    main()

