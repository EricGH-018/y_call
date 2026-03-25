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
from pathlib import Path

## IMPORT FUNCTIONS ##

from pulldown import load_snp_file_OY

def main():

    ## Load all markers in the Y chromosome ##
    OY_copy, OY, filtered_pos_OY = load_snp_file_OY() # import all SNPs in OY .csv file. 
    OY = OY.sort_values(by="pos") # sort by position.

    savepath = "./data/output/OY_snps.bed" # save a BED file with coordinates for every SNP. 

    dft = OY[["chrom", "pos"]].copy()
    dft["pos1"] = dft["pos"]
    dft.to_csv(savepath, sep="\t", index=False, header=None)
    print(f"Saved {len(dft)} OY SNPs to {savepath}")

    OY_copy.to_csv("data/output/merged_hg38_hg19_OY_SNP.csv",header=True)
    print(f"Saved correspondency between genome positions as: data/output/merged_hg38_hg19_OY_SNP.csv")

    filtered_pos_OY.to_csv("data/output/filtered_pos_OY.csv")
    print(f"Saved filtered positions as: data/output/filtered_pos_OY.csv")

if __name__ == "__main__":

    """ Run function main() only when the program is called from the CLI"""

    main()

