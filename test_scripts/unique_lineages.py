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

from pulldown import unique_lineages

def main():
    with open ("data/output/paths.txt", "r") as file:
        items = file.readlines()
    
    data = []
    for i in items:
        data.append(i.strip())


    summary_df = pd.read_csv("/home/eric_garcia_hoyos/y_call/data/output/scores.csv",header=None)
    summary_df.columns = ["Sample","Most derived","Level","Prop.","Score"]

    unique = unique_lineages(summary_df, data)

    print(f"Total number of unique lineages in dataset: {len(unique)}.")

if __name__ == "__main__":

    """ Run function main() only when the program is called from the CLI"""

    main()


