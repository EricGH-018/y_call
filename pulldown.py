## IMPORT MODULES ##

import numpy as np
import os  # for saving to folder
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcdefaults() # default settings to avoid warnings related to font.
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout

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

def create_parent_dct(path_parents):

    """ Return a dictionary for child (keys) and parent (values). """

    df_chpar = pd.read_csv(path_parents, header=None, usecols=[0, 1], names=["child", "parent"]) # store content of child - parent haplogroups in a dataFrame.
    chpar = dict(zip(df_chpar["child"], df_chpar["parent"])) # transform previous dataFrame into a dictionary.
    
    return chpar

def call_par(string, chpar):

    """ Return the level at which a specific haplogroup is found, by using a recursive function
    that checks every time the father of a child. """

    parent = chpar.get(string) # get the parent
    if parent is None: # if the consulted child is not inside the dictionary.
        return 0
    return call_par(parent, chpar) + 1 # if it is, add 1 and call again the function.

def ancder_par(string, par_dict, chpar):

    """ Return the total number of ancestral, derived and uncovered SNPs for parental branches """

    parent = chpar.get(string) # get the parent.
    if parent is None: # if the consulted child is not inside the dictionary.
        return 0
    return ancder_par(parent, par_dict, chpar) + par_dict.get(parent, 0) # if it is, add the number of SNPs for the parent and call again the function.

def support(string, chpar, derived_map):

    """ Return the total number of derived parental branches (known as support) """

    parent = chpar.get(string) # get the parent
    count = derived_map.get(parent, 0) # get the proportion of derived respect to total SNPs in branch. 
    
    if parent is None: # if the consulted child is not inside the dictionary.
        return count # return the count for the branch.
    
    return count + support(parent, chpar, derived_map) # if it is, increase the count by calling again the function.

def load_snp_file_OY(path_snps, reference_genome, chpar, branches, translation, unique=True):

    """ Return a dataFrame in Eigenstrat Format for SNPs,
    filtered for available positions, biallelic SNPs, Ref and Alt different,
    ACTG nucleotides and unique positions. Also, add information for level, 
    Y-haplogroup and YFull translation for every SNP. """

    # Read the content of the .csv file and assign column names.
    df_raw = pd.read_csv(path_snps, dtype=object, sep=" ", header=None, 
                         names=['Coordinates', 'ANC', 'DER', 'SNP-ID'])

    print(f"Loaded {len(df_raw)} SNPs")

    if reference_genome == "hg38":
        # Convert genome coordinates from hg38 to hg19
        chrY_len = 59373566 # chromosome length (in hg19)

        converter = get_lifter('hg38', 'hg19', one_based=True) # function to do the liftover.
        chrom = 'chrY'

        # Helper function to do the liftover (used then in map()).
        def convert_pos(p):
            conversion = converter[chrom][int(p)] # apply the liftover.
            return int(conversion[0][1]) if conversion else None # only return the position (if available).

        df_raw["pos"] = df_raw["Coordinates"].map(convert_pos) # add a new column named pos in df_raw, using map method to apply function convert_pos to every coordinate.
        df_raw = df_raw.drop(columns=['Coordinates']) # exclude column named Coordinates.

        # Filter for those positions out of the range for chrY (version hg19).
        filtered_df = df_raw["pos"] > chrY_len
        df_raw = df_raw[df_raw["pos"] < chrY_len] # define again df_raw but excluding rows where pos is out of range.
        print(f"# Out of range positions in hg19: {filtered_df.sum()}")

    elif reference_genome == "hg19":
        df_raw = df_raw.rename(columns={'Coordinates': 'pos'}) # if already in hg19 version, just rename the column Coordinates. 

    # Filter for non available position (NaN values). 
    idx = ~df_raw["pos"].isna()
    print(f"# Positions available: {np.sum(idx)}")

    df = df_raw[idx].reset_index(drop=True)
    df["pos"]=df["pos"].astype("int")

    # Filter only for biallelic - single based SNPs. 
    idx_bi= (df["ANC"].str.len()==1) & (df["DER"].str.len()==1) # one base on each column.
    print(f"# Biallelic SNPs: {np.sum(idx_bi)}")

    df = df[idx_bi].reset_index(drop=True)
    
    # Rename columns in dataFrame (useful in posterior analyses).
    df["ref"] = df["ANC"]
    df["alt"] = df["DER"]

    # Create a new column to indicate the chromosome ("Y"). 
    df["chrom"] = "Y"

    # Define the order for column names to be stored in the final df.
    cols = ["SNP-ID", "chrom", "pos", "ref", "alt"]
    df = df[cols]
    df = df.replace(regex=[' ','\n'], value='_')

    # Order SNPs by position.
    df = df.sort_values(by="pos")

    # Filter for SNPs with different alternative and derived states. 
    df = df[df["ref"] != df["alt"]]
    print(f"# Ref & Alt different: {len(df)}")

    # Filter for SNPs with only [ACTG] alleles. 
    snps_acceptable = ["A", "C", "T", "G"]
    df = df[df["ref"].isin(snps_acceptable) & df["alt"].isin(snps_acceptable)]
    print(f"# Ref & Alt ACTG: {len(df)}")

    # Filter for unique positions (non repeated coordinates). 
    if unique:
    	df = df.drop_duplicates(subset=["pos", "ref", "alt"], keep="first")
    	print(f"# Unique SNP positions: {len(df)}")

    print("Loading haplogroups information...")
    # Create a first dictionary to store values on haplogroups. 
    OY_dict = {}

    with open(branches, "r") as f:
        for line in f:
            items = line.strip().split() # because every line is: X-XXXX X123 X542 XG56*
            if not items:
                continue

            hap = items[0] # only the first value (the haplogroup).
            snps = items[1:] # the rest of values (the SNPs). 

            for snp in snps:
                if snp not in OY_dict:
                    OY_dict[snp] = []

                OY_dict[snp].append(hap) # every snp (key) with the corresponding haplogroup (value, repeated for some snps).

    # Create a second dictionary for Yfull translations in the identified haplogroups. 
    trans_dict = {}

    if pathlib.Path(translation).exists():
        with open(translation, "r") as f:
            for line in f:
                # split by whitespace
                items = line.strip().split(sep=",") # because every line is: A-XXXX,A-VWXXX
                if not items:
                    continue

                hapl = items[1]
                yfull_conv = items[0]

                trans_dict[hapl] = yfull_conv # every haplogroup (key) with its translation (value).


    # Once dictionaries are creted, consult every SNP and Haplogroup. 
    df["temp_hap"] = None
    df["temp_trans"] = None

    for index, row in df.iterrows():
        ids = str(row["SNP-ID"]).split(",") # because a SNP can have more than one ID: XUYDB,ZHJD,YODN.

        # Keep only unique entries as:
        found_haps = []
        for snp_id in ids: # check every name.
            h = OY_dict.get(snp_id) or OY_dict.get(snp_id+"*") # get the haplogroup for the SNP.
            if h and h not in found_haps: # if found and not already stored.
                found_haps.append(h)

        # If no haplogroup found, put None for haplogroup and its translation.
        if not found_haps:
            current_hap = ["None"]
            current_trans = ["None"]

        # If found:
        else:
            first_match = found_haps[0] # extract first element from the list, as there is only one).
            current_hap = first_match # get the series of haplogroups found (maybe only one). 
            current_trans = [trans_dict.get(h, "None") for h in current_hap] # get translation to YFull for every haplogroup.

        # Store lists in orignal dataFrame
        df.at[index, "temp_hap"] = current_hap
        df.at[index, "temp_trans"] = current_trans

    # Expand the dataFrame according to the lists stored in temporary columns
    df = df.explode(["temp_hap", "temp_trans"]).reset_index(drop=True)

    # Rename columns and delete temporary columns
    df["Y-haplogroup"] = df["temp_hap"] 
    df["YFull translation"] = df["temp_trans"]
    
    df = df.drop(columns=["temp_hap", "temp_trans"])

    # Make a search for the level of every haplogroup in the tree and add it as a new column.
    df["Level"] = df["Y-haplogroup"].map(lambda x: call_par(x, chpar))

    return df.reset_index(drop=True)

def ref_alt_count(df_ch, bases=["A", "C", "G", "T"]):

    """ Count Ref and Alt alleles in dataFrame df_ch
    with ref, alt, A, C, G, T fields and enter new columns
    ref# and alt# """

    # Define two new columns filled with 0.
    df_ch["ref#"]=0
    df_ch["alt#"]=0

    for base in bases:
        # Update ref# where the 'ref' column matches the current base
        df_ch.loc[df_ch["ref"] == base, "ref#"] = df_ch[base]
        
        # Update alt# where the 'alt' column matches the current base
        df_ch.loc[df_ch["alt"] == base, "alt#"] = df_ch[base]

    return df_ch

def pulldown_bamtable(path_bam = "", o_file = "",
                      bamtable = "./bin/BamTable",
                      snip5=0, snip3=0, base_qual=20, map_qual=25, 
                      path_bed = ""):

    """ Command for pileup using a .bam file and the corresponding coordinates in a .bed file, to create an output file. """

    run_cmd = f"{bamtable} -F -A --snip5={snip5} --snip3={snip3} --base_qual={base_qual} --map_qual={map_qual} -f {path_bed} {path_bam} > {o_file}"
    os.system(run_cmd)

def call_y_bam(df=[], path_bam="",
               path_bed = "",
               path_temp="",
               snip5=0, snip3=0, base_qual=20, map_qual=25):

    """ Return the pileup table after running the commad in a .bam file, both for all SNPs and derived sites only."""

    # Perform sanity checks whether input is there
    assert(os.path.exists(path_bam))
    assert(os.path.exists(path_bed))

    # Create the Pulldown
    pulldown_bamtable(path_bam = path_bam, o_file = path_temp, 
                      snip5=snip5, snip3=snip3, base_qual=base_qual, map_qual=map_qual,
                      path_bed = path_bed)

    # Load pulldown to read information contained in it.
    df1 = pd.read_csv(path_temp, sep="\t", header=None)
    df1.columns = ["chrom", "pos", "A", "C", "G", "T"] # read columns.

    # Select only columns for chrY and change it to only "Y". 
    idx = df1["chrom"]=="chrY"
    if np.sum(idx)>0:
        print(f"Changing {np.sum(idx)} ChrY -> Y")
        df1.loc[idx, "chrom"] = "Y"

    # Include information for the pulldown in the dataFrame passed to the function (comming from load_snp_file_OY())
    df2 = pd.merge(df, df1, on=["chrom", "pos"])

    # Show a series of coverage Statistics
    cov = df1[["A", "C", "G", "T"]].values
    cov1 = np.sum(cov, axis=1) # get the total count of reads per SNP.
    print(f"Average Coverage: {np.sum(cov1)/len(df):.4f}x") # sum the total of reads analysed and divide by the whole number of sites = Coverage as a proportion of the total.
    print(f"Sites covered: {np.sum(cov1>0)} / {len(df)}") # sites covered as those with a count higher than 0. 

    # Count the total number of reference and alternative alleles in each SNP. 
    df_ch = ref_alt_count(df2)

    # Get only SNPs with derived state (number of derived alleles higher than ancestral).
    idx_der = df_ch["alt#"]>df_ch["ref#"]
    print(f"Derived Loci: {np.sum(idx_der)} / {np.sum(cov1>0)} covered>0")

    # Create a dataFrame to store only derived SNPs.
    df_der = df_ch[idx_der].sort_values(by="SNP-ID").reset_index(drop=True).copy()

    return df_ch, df_der


def div_anc_der(df, chpar, df_exclude=[]):

    """ Return a dataFrame with total, ancestral, derived, and uncovered SNPs 
    per branch, as well as for parental branches. df_exclude: dataFrame of SNPs to filter """

    # If more than 1 SNP to be excluded in array, apply function to filter out the corresponding positions.
    if len(df_exclude)>0:
        df = exclude_snps(df, df_exclude)

    # Create columns for ANC and DER SNPs.
    df = df.copy()
    df["ancestral"] = df["alt#"] < df["ref#"]
    df["derived"]   = df["alt#"] > df["ref#"]

    # Divide the df by the different haplogroups.
    grouped = df.groupby("Y-haplogroup")

    # Generate a new df with columns for the Branch, the Level and the total for SNPs.
    new_df = grouped.agg(
        Level=("Level", "first"),
        Total_SNPs=("Y-haplogroup", "size"),
        Ancestral=("ancestral", "sum"),
        Derived=("derived", "sum"),
    ).reset_index()

    # Compute uncovered as the difference of the total with ancestral and derived.
    new_df["Uncovered"] = (
        new_df["Total_SNPs"]
        - new_df["Ancestral"]
        - new_df["Derived"]
    )

    # Order values in df by the score.
    new_df = new_df.sort_values(by="Level")

    # Use previously defined function to look for all ANC, DER and UNC SNPs in parental branches. First, transform into a dictionary to speed up processing.
    anc_lookup = new_df.set_index("Y-haplogroup")["Ancestral"].to_dict()
    der_lookup = new_df.set_index("Y-haplogroup")["Derived"].to_dict()
    unc_lookup = new_df.set_index("Y-haplogroup")["Uncovered"].to_dict()

    anc_par = []
    der_par = []
    unc_par = []

    # For each branch, find the total number of ancestral and derived SNPs in parent using ancder_par() function.
    for branch in new_df["Y-haplogroup"]:
        anc_par.append(ancder_par(branch, anc_lookup, chpar=chpar))
        der_par.append(ancder_par(branch, der_lookup, chpar=chpar))
        unc_par.append(ancder_par(branch, unc_lookup, chpar=chpar))

    # Insert new columns.
    new_df["#ANC in par."] = anc_par
    new_df["#DER in par."] = der_par
    new_df["#UNC in par."] = unc_par

    # Insert a new column as the total count of SNPs for parental branches.
    new_df["Total in par."] = new_df["#ANC in par."] + new_df["#DER in par."] + new_df["#UNC in par."]

    return new_df

def exclude_snps(df_ch, df_ex, col="SNP-ID"):

    """ Return dataFrame df_ch with filtered SNPs in df_ex. """

    idx= df_ch[col].isin(df_ex[col]) # look for specific SNP names.
    df_ch2 = df_ch[~idx].copy() # exclude rows for filtered SNPs. 

    return df_ch2

def create_tree(string, level, df, chpar, file=None, total_par=0,RED = "\033[91m", GREEN = "\033[92m", GREY = "\033[38;5;245m", RESET = "\033[0m"):

    """ Create a tree for an haplogroup in reverse order, and write it down in the corresponding file. """
    
    # Extract row for the haplotype considered. 
    subset = df.loc[df["Y-haplogroup"] == string]
    
    # If the haplogroup is found, extract the total number of ancestra and derived states, also in parental branches.
    if not subset.empty:

        row = subset.iloc[0]

        der = f"{GREEN}{row['Derived']}{RESET}"
        anc = f"{RED}{row['Ancestral']}{RESET}"
        der_par = f"{GREEN}{row['#DER in par.']}{RESET}"
        anc_par = f"{RED}{row['#ANC in par.']}{RESET}"
        total = f"{GREY}/{row['Total_SNPs']}{RESET}"

        total_par = int(row["Total in par."])
        string_par = f"{GREY}/{total_par}{RESET}"

    # If not found, add to the tree but with "Not found" strings.
    else:
        not_found = f"{GREY}0, Not Found{RESET}"
        der = anc = der_par = anc_par = not_found 
        string_par = total = ""

    # If a parent exists, use it to make a new call
    parent = chpar.get(string)
    if parent:
        create_tree(parent, level - 1, df, chpar, file=file, total_par=total_par) # here resting 1 to the level refers to the reverse order.

    # Create the tree starting at the root node.
    indent = "  " * level
    output = (
        f"{indent}|\n"
        f"{indent}|___> {string}, Level: {level}, "
        f"DER in branch: {der}{total}, "
        f"ANC in branch: {anc}{total}, "
        f"DER in par.: {der_par}{string_par}, "
        f"ANC in par.: {anc_par}{string_par}\n"
    )
    file.write(output)

def create_path(string, df, chpar, path=None):

    """ Create a path as [sample: XXX1, XXXX2, XX32...] """

    if path is None:
        path = []

    parent = chpar.get(string)

    if parent:
        if parent not in path:
            path = create_path(parent, df, chpar) # recursive step.

    path.append(string) # add node to the path.
    return(path)

def common_hap(paths):

    """ Return the common haplogroup between a subste of paths (wether it is the MRCA or the common
    most derived Y-haplogroup) """

    path = paths[0] # define path as the first entry of the list (starting point).

    for cons_path in paths[1:]: # for every path contained in the list provided to the function

        path = [x for x in path if x in cons_path] # extract only common elements between paths.

    # Define the haplogroup as the last entry of the path. 
    haplogroup = path[-1]

    return haplogroup

def network(data, width, height, avg_snps):

    """ Create a graphical network connecting Y-haplogroups and placing samples depending
    on the most derived Y-haplogroup assigned """

    G = nx.DiGraph()

    samples = set()

    for line in data: # here data as the file containing every of the paths to each sample.
        if not line or ":" not in line:
            continue
        
        # Store information for sample id and the different Y-haplogroups found.
        sample_id, path_str = line.split(':')
        samples.add(sample_id)
        nodes = [n.strip() for n in path_str.split(',') if n.strip() and n.strip().lower() != 'nan']
        
        # Create connections between nodes, with distance as the number of average_snps in the branch.
        for i in range(len(nodes) - 1):
            parent = nodes[i]
            child = nodes[i+1]

            snp_len = avg_snps.get(child, 1)

            G.add_edge(parent, child, minlen=snp_len)
        
        # Link sample to its specific leaf
        G.add_edge(nodes[-1], sample_id, len=1)
    
    # Force vertical layout
    G.graph['graph'] = {'rankdir':'TB'}
    pos = nx.nx_pydot.graphviz_layout(G, prog='dot')
    
    # Specify dimensions of the network
    plt.figure(figsize=(width, height)) 
    
    # Make a distinction between sample and haplogroup nodes.
    haplos = [n for n in G.nodes if n not in samples]
    
    # Draw background edges
    nx.draw_networkx_edges(G, pos, arrows=True, edge_color='black', alpha=0.2, width=0.5)
    
    # Create a dictionary of labels for edges that exceed the threshold, just to take into account the length of every branch.
    edge_labels = {}
    for u, v, data in G.edges(data=True):
        weight = data.get('minlen', data.get('len', 0))
        
        edge_labels[(u, v)] = f"{weight}"

    # Incorporate the branch length as a label in the network. 
    nx.draw_networkx_edge_labels(
        G, 
        pos, 
        edge_labels=edge_labels, 
        font_size=7,
        font_color='black',
        rotate=True,
        alpha=0.8,
        label_pos=0.5,
        verticalalignment = 'bottom',
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=0.5)
    )

    # Define also at the bottom of the network which samples correspond to every macrohaplogroup.
    roots_OY =  ["E-M35", "G-P15", "I-M253", "I-S238", "I-L416", "J-M267", "J-PAGES00028","L-M22", "R-L146", "R-M343", "T"]
    roots_ytree =  ["E-M35", "G-P15", "I1", "I-L460", "I-Y283553", "J1", "J2","L-M22", "R1a", "R1b", "T"]    
    roots = list(set(roots_OY + roots_ytree))

    # Set of colors for every macrohaplogroup.
    colors = {
        "E-M35": "#1b9e77",
        "G-P15": "#d95f02",
        "I-M253": "#7570b3",
        "I1": "#7570b3",
        "I-S238": "#e7298a",
        "I-L460": "#e7298a",
        "I-L416": "#7570b3",
        "I-Y283553": "#7570b3",
        "J-M267": "#e6ab02",
        "J1": "#e6ab02",
        "J-PAGES00028": "#a6761d",
        "J2": "#a6761d",
        "L-M22": "#90EE90",
        "R-L146": "#fb9a99",
        "R1a": "#fb9a99",
        "R-M343": "#e31a1c",
        "R1b": "#e31a1c",
        "T": "#6a3d9a"
    }
    
    # Alternative names of Y-haplogroups.
    alt_names = {
        "E-M35": "E1b1b",
        "G-P15": "G2a",
        "I-M253": "I1",
        "I-S238": "I2a",
        "I-L460": "I2a",
        "I-Y283553": "I2b",
        "J-M267": "J1",
        "J-PAGES00028": "J2",
        "L-M22": "L-M22",
        "R-L146": "R1a",
        "R-M343": "R1b",
        "T": "T"
    }
    
    # Extract all descendants nodes for a given macrohaplogroup.
    for node in haplos:
        if node in roots:
            downstream_nodes = nx.descendants(G, node)
            downstream_nodes.add(node) # include the starting node itself
    
            highlight_edges = []
    
            for u, v in G.edges():
                if u in downstream_nodes and v in downstream_nodes:
                    highlight_edges.append((u, v))
            
            # Draw edges in network for this path in the color stored inside stpolten_colors.
            nx.draw_networkx_edges(G, pos, edgelist=highlight_edges, 
                                   arrows=True, edge_color=colors[node], width=2.0, alpha=0.9)
    
    
    # Draw nodes for samples and haplogroups with high-density settings (smalle node sizes and font sizes).
    nx.draw_networkx_nodes(G, pos, nodelist=haplos, node_shape='s', 
                           node_color='#E0E0E0', node_size=200, edgecolors='grey')
    
    nx.draw_networkx_nodes(G, pos, nodelist=samples, node_shape='s', 
                           node_color='powderblue', node_size=500, edgecolors='black')
    
    # Label with a tiny font and background 'halo' for legibility
    nx.draw_networkx_labels(G, pos, font_size=5, font_family='sans-serif',
                            bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=0.2))

    y_values = [p[1] for p in pos.values()] # get all y values in the network.
    floor_y = min(y_values) # minimum y value.
    
    for node in samples:
        if node in pos:
            x, _ = pos[node]
            pos[node] = (x, floor_y)

    # Create a little offset for the highest value in the lables for macrohaplogroups.
    lower_offset = (max(y_values) - floor_y) * 0.05
    bar_y = floor_y - lower_offset 
    
    layer_y = 0 # define a layer depending on the macrohaplogroup (avoids overlapping)"
    for root in roots:
        if root not in G.nodes:
            continue
            
        # Get EVERY node in the sub-tree (samples AND internal haplos).
        sub_tree_nodes = nx.descendants(G, root)
        sub_tree_nodes.add(root)
        
        # Get X coordinates for every single node in this macrohaplogroup.
        x_coords = [pos[n][0] for n in sub_tree_nodes]
        
        if not x_coords:
            continue
            
        # Calculate x-range and add a small visual buffer (+/- 10 units)
        x_min, x_max = min(x_coords), max(x_coords)
        difference = x_max - x_min
        
        # Assign a color for the label
        color = colors.get(root, "#808080")
        
        # Draw the horizontal bar
        plt.hlines(y=bar_y + layer_y, xmin = x_min, xmax = x_max, 
                   color=color, linewidth=9)
        
        # Label the bar (slightly below).
        plt.text((x_min + x_max) / 2, (bar_y + layer_y) - 80, alt_names.get(root, root), 
                 color="black", ha='center', va='top', 
                 fontsize=30, rotation=0)

        layer_y -= 300 # adjust layer
    
    # Add title
    plt.title(f"Haplogroup Phylogeny: {len(samples)} Samples", fontsize=40)
    plt.axis('off')
    
    # Save as PDF and PNG.
    plt.savefig("figures/large_haplo_tree.pdf", bbox_inches='tight')
    plt.savefig("figures/large_haplo_tree.jpg", bbox_inches="tight")

def unique_lineages(df, data):

    """ Return a list of samples that correspond to unique lineages (paths shared by different samples counted as 1). """ 
    
    # Save all sample names in a list.
    samples = df["Sample"]
    samples = samples.to_list()

    # Store samples in data as a dictionary, to make a faster and more efficient search.
    path_map = {}
    for line in data:
        if ":" not in line:
            continue
        sample_id, path_str = line.split(':')
        nodes = {n.strip() for n in path_str.split(',') if n.strip() and n.strip().lower() != 'nan'}
        path_map[sample_id] = nodes

    # Define a list of unique samples
    unique_samples = set(df["Sample"])

    for samp, most_derived in zip(df["Sample"], df["Most derived"]):
        for other_id, other_nodes in path_map.items():

            # If it's a different sample AND the haplogroup is in their path
            if samp != other_id and most_derived in other_nodes:

                # If this sample is still in the set, remove it
                unique_samples.discard(samp)

                break # move to the next 'samp' once the current one is invalidated

    return list(unique_samples)

def y_call(bam_list, initial, final, base_qual, map_qual, database, reference_genome, create_network, width, height, transitions, translation, ex_limit, ages):

    """ Main function to run the y_call software. Other functions from the script included here."""

    out_path = pathlib.Path("data/output")
    out_path.mkdir(parents=True, exist_ok=True)

    # Make sure that the directory to the temporary file is created.
    temp_path = pathlib.Path("./data/output/temp")
    temp_path.mkdir(parents=True, exist_ok=True)

    path_snps=f"./data/input/{database}/snps.csv"
    tree=f"./data/input/{database}/tree.csv"
    mm=f"./data/input/{database}/mm.tsv"
    branches=f"./data/input/{database}/haps.csv"

    # Create a dictionary to store relations between child (key) and parent (value) haplogroups.
    chpar = create_parent_dct(path_parents=tree)

    print("\n## Process input data ##\n--------------------")
    # Path to the set of markers considered in the analysis (after filtering). 
    file_OY = pathlib.Path(f"./data/output/all_snps_{database}.csv")

    if file_OY.exists():
        print(f"Dataset for filtered markers already saved: ./data/output/all_snps_{database}.csv")
        OY = pd.read_csv(file_OY)
        path_bed = pathlib.Path(f"./data/output/{database}_snps.bed")
        if path_bed.exists:
            print(f"Corresponding .bed file in {path_bed}")
        else:
            savepath = f"./data/output/{database}_snps.bed" # save a BED file with coordinates for every SNP. 

            dft = OY[["chrom", "pos"]].copy()
            dft["pos1"] = dft["pos"]
            dft.to_csv(savepath, sep="\t", index=False, header=None)
            print(f"Coresponding .bed file in {path_bed}")
    else:
        print("Dataset for filtered markers missing\nCreating file...")
        OY = load_snp_file_OY(path_snps, reference_genome, chpar, branches, translation, unique=True) # import all SNPs in .csv file. 
        OY.to_csv(f"data/output/all_snps_{database}.csv")
        print(f"Done\nSaved as ./data/output/all_snps_{database}.csv")

        OY = OY.sort_values(by="pos") # sort by position.

        savepath = f"./data/output/{database}_snps.bed" # save a BED file with coordinates for every SNP. 

        dft = OY[["chrom", "pos"]].copy()
        dft["pos1"] = dft["pos"]
        dft.to_csv(savepath, sep="\t", index=False, header=None)
        print(f"Coresponding .bed file in {savepath}")

    print("\n## Load individuals ##\n--------------------")
    dft = pd.read_csv(bam_list, sep="\t") # file with all routes to BAM files.

    n = dft[dft["bam"].notna()] # select only rows where BAM file is available. 
    males = n[n["sex"] == "m"] # select only male samples.

    bam_dict_males = dict(zip(males["iid"], males["bam"])) # create a dictionary as [SAMPLE: route .bam].

    if final == 0:
        routes = list(bam_dict_males.items())[initial:final+len(bam_dict_males)+1] # if no indexes provided, consider all samples. 
    else:
        routes = list(bam_dict_males.items())[initial:final] # if indexes provided, select only those samples.     

    print(f"Loaded a total of {len(routes)} male individuals.")

    # Store in a list content for the final dataFrame:
    samples = []
    branches = []
    levels = []
    ratios = []
    scores = []

    # Store in lists ages for every haplogroup:
    ages_exist = False
    if pathlib.Path(ages).exists():
        formed = []
        tmrca = []
        ages_exist = True

        # Read the file that stores information for age on every haplogroup.
        ages = pd.read_csv(ages)

    # Create a dictionary to store the total number of SNPs per haplogroup:
    dict_snps = {}

    for sample in routes:

        # Show content for the sample being analysed.
        iid = sample[0]
        bam = sample[1]

        if not pathlib.Path(bam).exists(): # only if .bam file was not found.
            print(f"Sample {sample} not consdiered, missing BAM file.")
            continue

        print(f"\n## Summary coverage statistics for sample {iid} ##\n--------------------")

        # Do the call, depending on if the user has specified a minimum base quality and mapping quality.
        df_ch, df_der = call_y_bam(df=OY, path_bam=bam,
                                   path_bed=f'./data/output/{database}_snps.bed',
                                   path_temp=f'./data/output/temp/{iid}_temp.tsv',
                                   snip5=0, snip3=0, base_qual=base_qual, map_qual=map_qual)

        # Post-process filtering (subset of SNPs to be exluded in data/input/mm_v3.tsv)
        df_cts = pd.read_csv(mm, sep="\t")
        df_ex3 = df_cts[df_cts["count"]>ex_limit] # SNPs that are ancestral in test sample derived chains.


        # Filter for possible deaminations.
        if transitions == "Y":
            idx_trans = df_ch[((df_ch["ref"]=="C") & (df_ch["alt"]=="T")) | ((df_ch["ref"]=="G") & (df_ch["alt"]=="A"))].index # filter for SNPs where ref. alleles are C or G, and the alternative alleles T or A (possible deaminations).
            df_ch = df_ch.drop(index=idx_trans)

        # Create a dataFrame to know the total number of SNPs per branch, and the corresponding ancestral and derived counts. 
        dft = div_anc_der(df_ch, chpar=chpar, df_exclude=df_ex3)

        # Add a new column to the dataset as the ratio of derived SNPs in respect to the Total SNPs.
        dft["Derived/Total"] = round(dft["Derived"] / dft["Total_SNPs"],6)

        # Add a new column to the dataset as the ratio of #ANC in par. in #DER in par.
        dft["#ANC in par./#DER in par."] = round(dft["#ANC in par."] / dft["#DER in par."],6)

        # Add a new column to the dataset for the number of supporting derived nodes in each haplogroup.
        derived_nodes = (dft.set_index("Y-haplogroup")["Derived/Total"] > 0.95).astype(float).to_dict()
        dft["Support"] = dft["Y-haplogroup"].map(lambda x: support(x, chpar, derived_nodes))

        # Add a new column to the dataset as the score for every haplogroup = #DER in par. + Derived - 3* (#ANC in par. + Ancestral) - Uncovered - #UNC in par.
        dft["Score"] = round((dft["#ANC in par."]*-3)+(dft["#DER in par."])+(dft["Ancestral"]*-3)+(dft["Derived"])+(dft["Uncovered"]*-2)+(dft["#UNC in par."]*-2) + (dft["Support"]),6)

        # Get the most derived haplogroup according to score and level classification:
        df = dft[dft["Y-haplogroup"]!="I-P38"].sort_values(by="Score") # exclude haplogroup I-P38, as it's got a different name.

        most_der = df.iloc[-1, -1] # score for the most derived haplogroup.
        result = df[df.iloc[:, -1] == most_der] # select rows with that score.
        
        # If more than one row with the highest score:
        if len(result) > 1:
            paths = []

            # Define paths for every haplogroup with the highest score
            for haplo in result["Y-haplogroup"]:
                paths.append(create_path(string=haplo,df=df,chpar=chpar))

            # Decide the MRCA from the paths provided
            haplogr = common_hap(paths)
            if haplogr in df["Y-haplogroup"]:
                result = df[df["Y-haplogroup"]==haplogr]
            else:
                result = OY[OY["Y-haplogroup"]==haplogr]

        # Create directories for SNP and Y-haplogroups information.
        file_path = pathlib.Path(f"./data/output/markers/{iid}")
        file_path.mkdir(parents=True, exist_ok=True)

        file_path2 = pathlib.Path(f"./data/output/haplogroups/{iid}")
        file_path2.mkdir(parents=True, exist_ok=True)

        print(f"\n## Process output files for sample {iid} ##\n--------------------")
        # Write a .tsv file to consult SNPs for every branch.
        df_ch.to_csv(f"data/output/markers/{iid}/snps_{iid}.tsv", sep="\t", index=False)
        print(f"Saved data for pileup as data/output/markers/{iid}/snps_{iid}.tsv")

        # Write a .tsv file to consult only derived SNPs for every branch.
        df_der.to_csv(f"data/output/markers/{iid}/snps_der_{iid}.tsv", sep="\t", index=False)
        print(f"Saved data for pileup (considering only derived sites) as data/output/markers/{iid}/snps_der_{iid}.tsv")

        # Write a .tsv file to consult information for every Y-haplogroup called. 
        df.to_csv(f"data/output/haplogroups/{iid}/haplogroup_{iid}.tsv", sep="\t", index=False)
        print(f"Saved data for Y-haplgroup calls as data/output/haplogroups/{iid}/haplogroup_{iid}.tsv")

        # Write entries in the dictionary to see the total of SNPs per haplogroup.
        for haplo, total_snps in zip(df["Y-haplogroup"], df["Total_SNPs"]):
            if haplo in dict_snps:
                dict_snps[haplo].append(int(total_snps))
            else:
                dict_snps[haplo] = [int(total_snps)]

        # Create directory for the tree files created in each sample.
        file_path3 = pathlib.Path(f"data/output/trees/{iid}")
        file_path3.mkdir(parents=True, exist_ok=True)

        # Create a tree to know where is the sample located in Y-haplogroups
        with open(f"data/output/trees/{iid}/tree_output_{iid}.txt", "w") as f:
            f.write(f"Tree created for sample {iid}, with most derived Y-haplogroup: {result["Y-haplogroup"].iloc[0]}\n")
            create_tree(string = result["Y-haplogroup"].iloc[0],level = int(result["Level"].iloc[0]), df=df, chpar=chpar, file=f)
            f.write("\n")

        # Create a file to store all paths (for every sample), useful to then create the final tree
        file_path = "data/output/paths.txt"
        existing_entries = set()

        if os.path.exists("data/output/paths.txt"):
            with open(file_path, "r") as f:
                existing_entries = {line.strip() for line in f if line.strip()}

        with open(file_path, "a") as file:
            path = create_path(string=result["Y-haplogroup"].iloc[0], df=df, chpar=chpar)
            content = f"{iid}: {','.join(map(str, path))}"
    
            if content not in existing_entries:
                file.write(content + "\n")
                existing_entries.add(content)
        print(f"Y-haplogroup simplified tree saved as data/output/trees/{iid}/tree_output_{iid}.txt")

        print("\n## Update global output files ##\n--------------------\nPath to most derived haplogroup included in data/output/paths.txt")

        # Keep info for summary
        samples.append(iid)
        branches.append(result["Y-haplogroup"].iloc[0])
        levels.append(int(result["Level"].iloc[0]))

        if "#ANC in par./#DER in par." in result.columns and "Score" in result.columns:
            ratios.append(float(result["#ANC in par./#DER in par."].iloc[0]))
            scores.append(int(result["Score"].iloc[0]))
        else:
            ratios.append(0.0)
            scores.append(0)

        # Keep infor for haplogroup ages
        if ages_exist:
            formed_id = ages.set_index("Y-haplogroup")["Formed"].get(result["Y-haplogroup"].iloc[0])
            formed.append(formed_id)
            tmrca_id = ages.set_index("Y-haplogroup")["TMRCA"].get(result["Y-haplogroup"].iloc[0])     
            tmrca.append(tmrca_id)

    # Create a final dataFrame to store all measures for every sample.
    summary_df = pd.DataFrame({
        "Sample": samples,
        "Most derived": branches,
        "Level": levels,
        "Ratio": ratios,
        "Score": scores
    })
    #summary_df = summary_df.drop_duplicates(subset=['Sample'], keep='last') # in case the algorith runs a sample twice.

    # Write a .csv file with the summary.
    summary_df.to_csv("data/output/scores.csv", mode="a", header=False, index=False)
    print(f"Saved Y-haplogroup calls for {len(samples)} individuals as data/output/scores.csv")

    # Create also a dataFRame to store information on age.
    if ages_exist:
        ages_df = pd.DataFrame({
            "Sample": samples,
            "Most derived": branches,
            "Formed": formed,
            "TMRCA": tmrca
        })

        # Write a .csv file for ages.
        ages_df.to_csv("data/output/hap_ages.csv", mode="a", header=False, index=False)
        print(f"Saved Y-haplogroup ages for {len(samples)} individuals as data/output/hap_ages.csv")

    # Write a .csv file with SNPs per branch.
    data_list = []
    for haplo, snps in dict_snps.items():
        for snp in snps:
            data_list.append({"Haplogroup": haplo, "Total_SNPs": snp})

    # Create dataFrame and save.
    df_snps = pd.DataFrame(data_list)
    df_snps.to_csv("data/output/snps_branch.csv", mode="a", header=False, index=False)

    # Compute the total number of unique lineages for the dataset.
    # Use info for different paths in each sample and the corresponding scores. 
    print("\n## Run analysis to know the total number of unique lineages in dataset ##\n--------------------")
    with open ("data/output/paths.txt", "r") as file:
        items = file.readlines()

    data = []
    for i in items:
        data.append(i.strip())

    unique = unique_lineages(summary_df, data)

    print(f"Total number of unique lineages in dataset: {len(unique)}.")

    # Create a network to place samples in a Y-haplogroup tree (only if specified by the user).
    if create_network == "Y":
        print("\n## Create network to know relationship between samples ##\n--------------------")
        df = pd.read_csv("data/output/snps_branch.csv", header=None)

        # Divide the dataFrame by the different haplogroups.
        df.columns = ["Haplogroup", "SNPs_num"]
        grouped = df.groupby("Haplogroup")

        # Generate a new dataFarme with columns for the Branch, the Level and the average total of SNPs.
        new_df = grouped.agg(
                Sum=("SNPs_num", "mean")
        ).round(2).reset_index()

        avg_snps = dict(zip(new_df["Haplogroup"], new_df["Sum"]))

        # Create directories to store figures.
        file_fig = pathlib.Path(f"./figures")
        file_fig.mkdir(parents=True, exist_ok=True)

        network(data, width, height, avg_snps)

        print(f"Tree created for a total of {len(data)} samples")

