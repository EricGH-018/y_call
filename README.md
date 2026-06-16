# y_call
A Python package designed to make Y-haplogroup calls using large sets of SNPs and adressed to deal with ancient DNA issues.

## Overview
`y_call` is a reproducible pipeline for studying Y‑haplogroups inferred from low‑coverage or ancient DNA BAM files, 
since it evaluates a large number of genomic positions in the Non-recombining Region of the Y chromosome (NRY). 
It performs SNP filtering, genotyping via pileup extraction, deamination filtering, and Y-haplogroup calls, 
placing samples into known Y-phylogenies and displaying a phylogenetic tree if desired by the user.

## Features

- SNP filtering and database‑driven Y-haplogroup structure.
- Pileup extraction and genotyping using BamTable (v0.1.2).
- Ancestral vs. derived allele quantification.
- Y-haplogroup scoring and MRCA resolution.
- Per‑sample tree generation.
- Global phylogenetic tree (optional).
- Summary statistics and lineage analysis.
- CLI and Python API.

## Installation

Install from PyPI:

```bash
pip install y_call
```

Also, if preferred, install from source:

```bash
git clone https://github.com/EricGH-018/y_call
cd y_call
pip install .
```

## Downloading reference databases

The `y_call` package does not bundle the full Y-chromosome reference datasets (YFull and OY) 
to keep the installation lightweight.  
Instead, these databases are hosted as versioned assets on GitHub Releases.

You can download them directly using the built‑in command:

```bash
y_call_download_data --database YFull
y_call_download_data --database OY
y_call_download_data --database translation
y_call_download_data --database ages
```
The corresponding files will be stored in:
```bash
~/data/input/
```
NOTE: always make sure that the working directory includes the `./data/input` directory, 
since `y_call` needs it to look for input data.

The YFull and OY reference datasets included in this package are derived from
publicly available information provided by YFull/YTeam and OY, respectively. 
These files are reformatted and processed versions created specifically for use with `y_call`.

Original data sources:
- YFull YTree: https://github.com/YFullTeam/YTree
- OY: https://ysnp.info/downloads/

The `y_call` package does not redistribute the original tree files.
Only derived, reformatted tables are included.


## Updating the file mm.tsv for filtering

`y_call` makes a filtering of certain SNPs found to exhibit ancestral states 
in 1 sample or more. 
After a first Y-haplofgroup call for a set of individuals, the user can update 
the corresponding file with:

```bash
y_call_filter_mm --database YFull
y_call_filter_mm --database OY
```
The existing mm.tsv considers this filtering from a subset of 247 samples.

## Unifying all Y-haplogroup calls in summary files

After Y-haplogroup calling at the sample level, `y_call` enables the user to create
file integrating information for the whole dataset analysed. This can be done with:
```bash
y_call_unify --database YFull --create_phylogeny Y
y_call_unify --database OY --create_phylogeny Y
``` 

## Command‑line usage (with default parameters)

```bash
y-call \
  --bam-list samples.tsv \
  --index_i 0 \
  --index_f 0 \
  --base_qual 20 \
  --map_qual 25 \
  --rg "hg38" \
  --width 60 \
  --height 80 \
  --transitions "N" \
  --k 3 \
  --database "YFull" \
  --translation "./data/input/YF-translations.csv"
  --ex_limit 5 \
  --ages "./data/input/ages.csv" \
  --snip3 0 \
  --snip5 0
```

## Python API (with default parameters)

```python
from pulldown import y_call

y_call(
    bam_list="samples.tsv",
    initial=0,
    final=0,
    base_qual=20,
    map_qual=25,
    k=3,
    database="YFull",
    reference_genome="hg38",
    width=60,
    height=80,
    transitions="N",
    translation="./data/input/YF-translations.csv",
    ex_limit=5,
    ages="./data/input/ages.csv",
    snip3=0,
    snip5=0
)
```

## Input data

From the root directory, `y_call` expects the following organisation:

```text
data/
  input/
    <database>/
      snps.csv
      tree.csv
      mm.tsv
      haps.csv
    ages.csv
    YF-translations.csv
```

## Output data

After execution, `y_call` generates the following output structure:
```text
data/
    output/
        markers/<idd>/
            snps_<iid>.tsv
            snps_der_<iid>.tsv
        haplogroups/<iid>/
            haplogroup_<iid>.tsv
        trees/<iid>/
            tree_output_<iid>.txt
        temp/
            <iid>.tsv
        hap_ages.csv
        paths.txt
        scores.csv
        simple_tree.nwk
        snps_branch.csv
figures/
    Y-haplogroup_tree.jpg
    Y-haplogroup_tree.pdf
```

To better understand the content of every output file, refer to the Jupyter Vignette 
stored in: https://github.com/EricGH-018/y_call/blob/master/Jupyter_Vignette/y_call_vignette.ipynb. 
The different files genertaed contain:

- `data/output/markers/<iid>/snps_<iid>.tsv` — SNP information: SNP-ID, chrom, pos, ref allele, 
alt allele, Y-haplogroup, YFull translation, Level, ACGT read counts, ref and alt read counts.
- `data/output/markers/<iid>/snps_<iid>.tsv` — SNP information (only for variants where 
alt count > ref count): SNP-ID, chrom, pos, ref allele, alt allele, Y-haplogroup, YFull translation, 
Level, ACGT read counts, ref and alt read counts.
- `data/output/haplogroups/<iid>/haplogroup_<iid>.tsv` — haplogroup scoring table: Y-haplogroup, 
Level, Total number of SNPs, SNPs supporting ancestral state, SNPs supporting alternative state, 
SNPs uncovered, SNPs supporting ancestral alleles in nodes until root, SNPs supporting alternative 
alleles in nodes until root, SNPs uncovered until root, Total number of SNPs until root, 
Proportion of derived SNPs compared to Total of SNPs, Ratio of ancestral SNPs and derived SNPs, 
Supporting nodes until Y-haplogroup and Score. 
- `data/output/trees/<iid>/tree_output_<iid>.txt` — tree file for an individual sample, using 
ASCII characters and indicating Alternative, Ancestral, Total SNPs, Alternative in par., 
Ancestral in par., Total in par. or nodes Not Found.
- `data/output/hap_ages.csv` — Formed and TMRCA ages for every Y-haplogroup call.
- `data/output/paths.csv` — List of nodes until root for assigned Y-haplogroups. 
- `data/output/scores.txt` — Summary statistics for Y-haplogroup calls: iid, Y-haplogroup, 
Level, Ratio Anc in par./Der in par., Score. 
- `data/output/simple_tree.nwk` — phylogenetic tree in newick format with branch lengths 
as average number of alternative SNPs positions.
- `data/output/snps_branch.csv` — total number of alternative SNPs found per node.
- `figures/Y-haplogroup_tree.*` — image of the phylogenetic tree (in .jpg and .pdf) in simple_tree.nwk.


