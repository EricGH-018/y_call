## IMPORT MODULES ##

import argparse
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import pandas as pd

## IMPORT FUNCTIONS ##

from pulldown import network

def main():

    parser = argparse.ArgumentParser(
        description = "Assign most derived Y haplogroups in samples."
    )

    parser.add_argument(
        "--width",
        required=False,
        help = "Width of the network created."
    )

    parser.add_argument(
        "--height",
        required=False,
        help = "Height of the network created."
    )


    # Store arguments as variables:
    args = parser.parse_args()

    width = args.width
    height = args.height

    with open ("data/output/paths.txt", "r") as file:
        items = file.readlines()
    
    data = []
    for i in items:
        data.append(i.strip())

    df = pd.read_csv("data/output/snps_branch.csv", header=None)
    
    # Divide the df by the different haplogroups.
    df.columns = ["Haplogroup", "SNPs_num"]
    grouped = df.groupby("Haplogroup")
    
    # Generate a new df with columns for the Branch, the Level and the total for SNPs.
    new_df = grouped.agg(
            Sum=("SNPs_num", "mean")
    ).round(2).reset_index()
    
    avg_snps = dict(zip(new_df["Haplogroup"], new_df["Sum"]))

    network(data, int(width), int(height), avg_snps)

if __name__ == "__main__":

    """ Run function main() only when the program is called from the CLI"""

    main()

