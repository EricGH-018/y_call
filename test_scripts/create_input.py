# IMPORT MODULES
import pandas as pd
import json
import re
import argparse

from pulldown import get_tree, get_haplog, get_age 

# Add a description for the script and the corresponding arguments passed to it.
parser = argparse.ArgumentParser(description = "Create necessary input files from a given tree in JSON format")

parser.add_argument("--database", required=True, type=str,
                    help = "Specify the dataset to be used in the analysis, between OY, YFull or other (the program looks for ./data/input/[database]/).")

parser.add_argument("--json", required=True, type=str,
                    help = "Give the path to a JSON file containing the tree.")

args = parser.parse_args()
database = args.database
json = args.json

with open(json, "r") as f:
    data = json.load(f)

tree_list = create_tree(data)
split_tree = [item.split(",") for item in tree_list]
df_tree = pd.DataFrame(split_tree, columns=["Parent", "Child"])
df_tree = df_tree[["Child", "Parent"]]
df_tree.to_csv(f"./data/input/{database}/tree.csv",index=False,header=False)

snps_list = get_haplog(data)
new_snps = []
for line in snps_list:
    line = re.sub(" ","",line)
    line = re.sub(","," ",line)
    line = re.sub("/"," ",line)
    new_snps.append(line)
with open("/home/eric_garcia_hoyos/y_call/data/input/haps.csv", "w") as file:
    for line in new_snps:
        file.write(line + "\n")

ages = get_age(data)
split_ages = [item.split(",") for item in ages]
df_age = pd.DataFrame(split_ages, columns=["Y-haplogroup","Formed", "Formed_low","Formed_high","TMRCA","TMRCA_low","TMRCA_high"])
df_age.to_csv("./data/input/ages.csv",index=False)
