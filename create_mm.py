# IMPORT MODULES
import argparse
import pandas as pd
from pulldown import create_parent_dct, get_mismatch_snps

# Add a description for the script and the corresponding arguments passed to it.
parser = argparse.ArgumentParser(description = "Update mm.tsv file if a great number of mismatches is found in Y-haplogroup calls.")


parser.add_argument("--database", required=False, type=str, default="YFull",
                                   help = "Specify the dataset to be used in the analysis, between OY, YFull or other (the program looks for ./data/input/[database]/).")
args = parser.parse_args()
database = args.database


print("\n--- RUNNING ANALYSIS TO UPDATE mm.tsv ---\n")

y_dct = {}

print("Reading file...")
chpar = create_parent_dct(path_parents=f"./data/input/{database}/tree.csv") # create a dictionary to know parent - child relations between haplogroups.
summary_df = pd.read_csv("./data/output/scores.csv", header=None, names = ["Sample_ID", "Y-haplogroup", "Level", "#ANC in par./#DER in par.", "Score"]) # get information from the Y-haplogroup calls.
print("Done")

# Create a new dictionary with every sample and its haplogroup: 
for _, row in summary_df.iterrows():
    y_dct[row["Sample_ID"]] = row["Y-haplogroup"]

df_mms_lst  = []

print("\nAnalysing information for every individual...")
for iid, hpgrp in y_dct.items():

    # For every individual, extract information for the pulldown:
    df_ch = pd.read_csv(f"./data/output/markers/{iid}/snps_{iid}.tsv", sep="\t")
    
    # Apply function to get those snps with a higher count of ancestral over derived.
    df_mms = get_mismatch_snps(hpgrp, chpar=chpar, df_ch=df_ch)

    # Store in list
    df_mms_lst.append(df_mms)

df_mms = pd.concat(df_mms_lst) # create a dataFrame from list.
df_cts = df_mms["SNP-ID"].value_counts() # only extract the times that a certain SNP has appeared.
df_cts.to_csv(f"./data/input/{database}/mm.tsv", sep="\t") # store in a .csv file.
print(f"Done\nSaved at ./data/input/{database}/mm.tsv")
