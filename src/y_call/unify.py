# IMPORT MODULES
import argparse
from html import parser
import pandas as pd
import os
import pathlib
from .pulldown import create_parent_dct, final_summary, final_ages, snps_branch, create_path, unique_lineages, generate_simple_grouped_tree, nwk_tree

def unify(database: str, create_phylogeny: str, width: int, height: int):

    # Empty lists to store the information for every sample: 
    samples = []
    branches = []
    levels = []
    ratios = []
    scores = []
    formed = []
    tmrca = []
    dict_snps = {}

    tree = f"./data/input/{database}/tree.csv"
    chpar = create_parent_dct(path_parents=tree)
    paths_file = pathlib.Path("./data/output/paths.txt")

    if pathlib.Path("./data/input/ages.csv").exists():
        ages_exist = True
        ages_df = pd.read_csv("./data/input/ages.csv")
    else:
        ages_df = ""

    base = "./data/output/haplogroups/"

    for name in os.listdir(base):
        full = os.path.join(base, name)
        sample = full.split("/")[-1]
        if os.path.isdir(full):
            haplo_call = pd.read_csv(f"{full}/haplogroup_{sample}.tsv", sep="\t")
            samples.append(sample)
            branches.append(haplo_call.tail(1)["Y-haplogroup"].iloc[0])
            levels.append(haplo_call.tail(1)["Level"].iloc[0])
            ratios.append(haplo_call.tail(1)["#ANC in par./#DER in par."].iloc[0])
            scores.append(haplo_call.tail(1)["Score"].iloc[0])

            for haplo, total_snps in zip(haplo_call["Y-haplogroup"], haplo_call["Derived"]):
                dict_snps.setdefault(haplo, []).append(int(total_snps))

            existing_entries = set()

            if paths_file.exists():
                with open(paths_file, "r") as f:
                    existing_entries = {line.strip() for line in f if line.strip()}

            path = create_path(
                string=haplo_call.tail(1)["Y-haplogroup"].iloc[0],
                df=haplo_call,
                chpar=chpar,
            )
            
            content = f"{sample}: {','.join(map(str, path))}"

            if content not in existing_entries:
                with open(paths_file, "a") as file:
                    file.write(content + "\n")            

            if ages_exist:
                ages_indexed = ages_df.set_index("Y-haplogroup")
                formed.append(ages_indexed["Formed"].get(haplo_call.tail(1)["Y-haplogroup"].iloc[0]))
                tmrca.append(ages_indexed["TMRCA"].get(haplo_call.tail(1)["Y-haplogroup"].iloc[0]))
    
    print(
            "\n## Update global output files ##\n--------------------\n"
            "Paths to most derived Y-haplogroups included in data/output/paths.txt"
        )
    final_summary(samples, branches, levels, ratios, scores)
    snps_branch(dict_snps)
    final_ages(ages_exist, samples, branches, formed, tmrca)

    if paths_file.exists():
        with open(paths_file, "r") as file:
            items = [line.strip() for line in file if line.strip()]
    else:
        items = []
    if create_phylogeny == "Y" and items:
        print(
            "\n## Create a phylogenetic tree to know relationship between samples ##\n"
            "--------------------"
        )
        avg_snps = "./data/output/snps_branch.csv"
        nwk_out = "./data/output/simple_tree.nwk"

        file_fig = pathlib.Path("./figures")
        file_fig.mkdir(parents=True, exist_ok=True)

        nwk_tree(paths_file, width, height, avg_snps, nwk_out, items)

        print(f"Tree created for a total of {len(items)} samples")

    print(
        "\n## Run analysis to know the total number of unique lineages in dataset ##\n"
        "--------------------"
    )

    if paths_file.exists():
        with open(paths_file, "r") as file:
            items = [line.strip() for line in file if line.strip()]
    else:
        items = []

    scores_file = pathlib.Path("./data/output/scores.csv")
    if scores_file.exists():
        summary = pd.read_csv(
            scores_file,
            header=None,
            names=[
                "Sample-ID",
                "Y-haplogroup",
                "Level",
                "#ANC in par./#DER in par.",
                "Score",
            ],
        )
        unique = unique_lineages(summary, items)
        print(f"Total number of unique lineages in dataset: {len(unique)}.")
    else:
        print("No scores.csv found; cannot compute unique lineages.")
        unique = []



def main():
    # Add a description for the script and the corresponding arguments passed to it.
    parser = argparse.ArgumentParser(description = "Create files to store global inforation on assigned Y-haplogroups and related data.")

    parser.add_argument("--database", required=True, type=str,
                                   help = "Select the database name to update the global files.")

    parser.add_argument("--create_phylogeny", required=True, type=str,
                                   help = "Create a network connecting the different Y-haplogroups and placing samples in it, specify Y if wanted.")

    parser.add_argument("--width", required=False, type=int, default=60,
        help = "Width of the network created.")

    parser.add_argument("--height", required=False, type=int, default=80,
        help = "Height of the network created.")

    args = parser.parse_args()

    unify(args.database, args.create_phylogeny, args.width, args.height)