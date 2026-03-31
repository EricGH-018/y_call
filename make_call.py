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

from pulldown import y_call, load_snp_file_OY, ref_alt_count, pulldown_bamtable, call_y_bam, create_tree, div_anc_der, create_parent_dct, create_path, network, unique_lineages

def main():

    # Define parameters for the algorithm.  
    parser = argparse.ArgumentParser(description = "Assign most derived Y haplogroups in mapping samples.")
    parser.add_argument("--bam_list", required=True, type=str,
        help = "Provide a .tsv file with valid routes to .bam files, columns: iid, bam, sex (f or m).")
    parser.add_argument("--index_i", required=False, type=int, default=0,
        help = "Specify which is the first sample in the set to analyse.")
    parser.add_argument("--index_f", required=False, type=int, default=0,
        help = "Specify which is the last sample in the set to analyse.")
    parser.add_argument("--base_qual", required=False, type=int, default=20,
        help = "Minimum base quality.")
    parser.add_argument("--map_qual", required=False, type=int, default=25,
        help = "Minimum mapping quality.")
    parser.add_argument("--snps", required=False, type=str, default="./data/input/YBrowse_snps_hg38-runtime.csv",
                        help = "Provide a .csv file containing information for SNPs in Y chromosome, including the following columns: pos, ref, alt, ID.")
    parser.add_argument("--rg", required=False, type=str, default="hg38",
        help = "Specify the reference genome for genome coordinates (either hg38 or hg19).")
    parser.add_argument("--create_network", required=False, type=str, default="N",
        help = "Create a network connecting the different Y-haplogroups and placing samples in it, specify Y if wanted.")
    parser.add_argument("--width", required=False, type=int, default=120,
        help = "Width of the network created.")
    parser.add_argument("--height", required=False, type=int, default=90,
        help = "Height of the network created.")
    parser.add_argument("--transitions", required=False, type=str, default="N",
        help = "Filter for transitions in the analysis (C -> T and G -> A), specify Y if wanted.")
    parser.add_argument("--haplogr_tree", required=False, type=str, default="./data/input/OYchpar.csv",
        help = "Provide a .csv file indicating for every haplogroup its parent. Two columns, first child node, and second parent node.")
    parser.add_argument("--translation", required=False, type=str, default="./data/input/YF-translations.csv",
        help = "Provide a .csv file indicating the equivalent nomenclature for every haplogroup in a different database.")
    parser.add_argument("--mismatches", required=False, type=str, default="./data/input/mm_v3.tsv",
        help = "Set of SNPs to be filtered because of high count of ancestral alleles in a great number of samples.")
    parser.add_argument("--branches", required=False, type=str, default="./data/input/OYhap.db",
        help = "Provide a .csv file with all SNPs contained in each haplogroup.")
    parser.add_argument("--ex_limit", required=False, type=int, default=5,
        help = "Specify the minimum number of samples to consdier excluding a SNP (because of high ancestral state count).")

    # Store arguments as variables:
    args = parser.parse_args()

    print("\n--- RUNNING y_call SOFTWARE ---\n")

    # Run function for the software: y_call.
    y_call(bam_list=args.bam_list, initial=args.index_i, final=args.index_f,
           base_qual=args.base_qual, map_qual=args.map_qual, path_snps=args.snps,
           reference_genome=args.rg, create_network=args.create_network, 
           width=args.width, height=args.height, transitions=args.transitions,
           tree=args.haplogr_tree, translation=args.translation, mm=args.mismatches, 
           branches=args.branches, ex_limit=args.ex_limit)

if __name__ == "__main__":

    """ Run function main() only when the program is called from the CLI"""

    main()

