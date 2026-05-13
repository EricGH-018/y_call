""" Python script storing test functions to check the correct functioning for the package """

import pandas as pd
import numpy as np
from unittest.mock import patch, mock_open, MagicMock
from y_call.pulldown import create_parent_dct, call_par, ancder_par, support,ref_alt_count,pulldown_bamtable, load_snp_file_OY, exclude_snps, div_anc_der, create_tree, create_path, common_hap, recurse, generate_simple_grouped_tree, unique_lineages,

@patch("y_call.pulldown.pd.read_csv")
def test_create_parent_dct_basic(mock_read_csv):
    
    # Create a fake .csv file
    df_fake = pd.DataFrame({
        "child": ["A", "B", "C"],
        "parent": ["X", np.nan, "Z"]
    })
    mock_read_csv.return_value = df_fake

    result = create_parent_dct("fake_path.csv")

    expected = {
        "A": "X",
        "B": "Y-Chromosome Adam",  # NaN replaced
        "C": "Z"
    }

    assert result == expected

def test_call_par_basic():
    chpar = {
        "A": "B",
        "B": "C",
        "C": "Y-Chromosome Adam"
    }
    assert call_par("A", chpar) == 3
    assert call_par("B", chpar) == 2
    assert call_par("C", chpar) == 1

def test_call_par_unknown():
    chpar = {"A": "B"}
    assert call_par("X", chpar) == 0

def test_ancder_par_basic():
    chpar = {
        "A": "B",
        "B": "C",
        "C": "Y-Chromosome Adam"
    }
    par_dict = {
        "B": 5,
        "C": 3,
        "Y-Chromosome Adam": 2
    }

    # A → B → C → Adam = 5 + 3 + 2 = 10
    assert ancder_par("A", par_dict, chpar) == 10

def test_ancder_par_unknown():
    chpar = {"A": "B"}
    par_dict = {"B": 5}
    assert ancder_par("X", par_dict, chpar) == 0

def test_support_basic():
    chpar = {
        "A": "B",
        "B": "C",
        "C": "Y-Chromosome Adam"
    }
    derived_map = {
        "B": 1,
        "C": 2,
        "Y-Chromosome Adam": 3
    }

    # A → B → C → Adam = 1 + 2 + 3 = 6
    assert support("A", chpar, derived_map) == 6

def test_support_unknown():
    chpar = {"A": "B"}
    derived_map = {"B": 1}
    assert support("X", chpar, derived_map) == 0

def test_ref_alt_count_basic():
    df = pd.DataFrame({
        "ref": ["A", "C"],
        "alt": ["G", "T"],
        "A": [5, 0],
        "C": [0, 3],
        "G": [2, 0],
        "T": [0, 7]
    })

    df_out = ref_alt_count(df.copy())

    assert df_out.loc[0, "ref#"] == 5
    assert df_out.loc[0, "alt#"] == 2
    assert df_out.loc[1, "ref#"] == 3
    assert df_out.loc[1, "alt#"] == 7

@patch("y_call.pulldown.os.system")
def test_pulldown_bamtable_builds_correct_command(mock_system):
    pulldown_bamtable(
        path_bam="sample.bam",
        o_file="out.txt",
        bamtable="./bin/BamTable",
        snip5=5,
        snip3=3,
        base_qual=20,
        map_qual=25,
        path_bed="coords.bed"
    )

    expected_cmd = (
        "./bin/BamTable -F -A --snip5=5 --snip3=3 "
        "--base_qual=20 --map_qual=25 -f coords.bed sample.bam > out.txt"
    )

    mock_system.assert_called_once_with(expected_cmd)

@patch("y_call.pulldown.pd.read_csv")
@patch("y_call.pulldown.pathlib.Path.exists")
@patch("y_call.pulldown.call_par")
def test_load_snp_file_basic(mock_call_par, mock_exists, mock_read_csv):
    # Fake CSV input
    df_fake = pd.DataFrame({
        "Coordinates": ["100", "200"],
        "ANC": ["A", "C"],
        "DER": ["G", "T"],
        "SNP-ID": ["rs1", "rs2"]
    })
    mock_read_csv.return_value = df_fake

    # Fake translation file exists
    mock_exists.return_value = False

    # Fake parent dictionary
    chpar = {"A": "B", "B": "C", "C": "Y-Chromosome Adam"}

    # Fake branches file
    fake_branches = "rs1 A\nrs2 B\n"
    with patch("builtins.open", mock_open(read_data=fake_branches)):
        df_out = load_snp_file_OY(
            path_snps="fake_snps.txt",
            reference_genome="hg19",
            chpar=chpar,
            branches="fake_branches.txt",
            translation="fake_translation.txt"
        )

    # Basic structure checks
    assert list(df_out.columns) == [
        "SNP-ID", "chrom", "pos", "ref", "alt",
        "Y-haplogroup", "YFull translation", "Level"
    ]

    assert len(df_out) == 2
    assert df_out.loc[0, "chrom"] == "Y"
    assert df_out.loc[0, "ref"] == "A"
    assert df_out.loc[0, "alt"] == "G"

@patch("y_call.pulldown.pd.read_csv")
@patch("y_call.pulldown.pathlib.Path.exists")
@patch("y_call.pulldown.call_par")
def test_load_snp_file_haplogroups(mock_call_par, mock_exists, mock_read_csv):
    df_fake = pd.DataFrame({
        "Coordinates": ["100"],
        "ANC": ["A"],
        "DER": ["G"],
        "SNP-ID": ["rsX"]
    })
    mock_read_csv.return_value = df_fake
    mock_exists.return_value = False

    fake_branches = "A rsX\n"
    with patch("builtins.open", mock_open(read_data=fake_branches)):
        df_out = load_snp_file_OY(
            path_snps="fake",
            reference_genome="hg19",
            chpar={"A": "B"},
            branches="fake_branches.txt",
            translation="fake_translation.txt"
        )

    assert df_out.loc[0, "Y-haplogroup"] == "A"

@patch("y_call.pulldown.pd.read_csv")
@patch("y_call.pulldown.pathlib.Path.exists")
@patch("y_call.pulldown.call_par")
def test_load_snp_file_translation_and_level(mock_call_par, mock_exists, mock_read_csv):
    df_fake = pd.DataFrame({
        "Coordinates": ["100"],
        "ANC": ["A"],
        "DER": ["G"],
        "SNP-ID": ["rsX"]
    })
    mock_read_csv.return_value = df_fake

    mock_exists.return_value = True  # translation file exists

    fake_branches = "A rsX\n"
    fake_translation = "A-V1,A\n"

    mock_call_par.return_value = 3

    with patch("builtins.open", mock_open(read_data=fake_branches)) as m1:
        with patch("builtins.open", mock_open(read_data=fake_translation)) as m2:
            df_out = load_snp_file_OY(
                path_snps="fake",
                reference_genome="hg19",
                chpar={"A": "B"},
                branches="fake_branches.txt",
                translation="fake_translation.txt"
            )

    assert df_out.loc[0, "YFull translation"] == "A-V1"
    assert df_out.loc[0, "Level"] == 3

def test_call_y_bam_pipeline(mock_exists, mock_pulldown, mock_read_csv, mock_ref_alt):
    mock_exists.return_value = True

    df_input = pd.DataFrame({"chrom": ["Y"], "pos": [100], "SNP-ID": ["rs1"]})

    df_pulldown = pd.DataFrame({
        "chrom": ["chrY"],
        "pos": [100],
        "A": [5], "C": [0], "G": [0], "T": [0]
    })
    mock_read_csv.return_value = df_pulldown

    df_ch_fake = pd.DataFrame({
        "chrom": ["Y"],
        "pos": [100],
        "SNP-ID": ["rs1"],
        "ref#": [5],
        "alt#": [0]
    })
    mock_ref_alt.return_value = df_ch_fake

    df_ch, df_der = call_y_bam(
        df=df_input,
        path_bam="fake.bam",
        path_bed="fake.bed",
        path_temp="temp.txt"
    )

    mock_pulldown.assert_called_once()
    assert len(df_ch) == 1
    assert len(df_der) == 0

def test_exclude_snps_basic():
    df_ch = pd.DataFrame({"SNP-ID": ["rs1", "rs2", "rs3"]})
    df_ex = pd.DataFrame({"SNP-ID": ["rs2"]})

    df_out = exclude_snps(df_ch, df_ex)

    assert list(df_out["SNP-ID"]) == ["rs1", "rs3"]

@patch("y_call.pulldown.ancder_par")
def test_div_anc_der_basic(mock_ancder):
    mock_ancder.return_value = 10  # constant for simplicity

    df_ch = pd.DataFrame({
        "Y-haplogroup": ["A", "A", "B"],
        "Level": [1, 1, 2],
        "ref#": [5, 0, 3],
        "alt#": [0, 2, 1]
    })

    df_ex = pd.DataFrame({"SNP-ID": []})  # no exclusions

    df_out = div_anc_der(df_ch, chpar={"A": "B", "B": None}, df_exclude=df_ex)

    # Check structure
    assert list(df_out.columns) == [
        "Y-haplogroup", "Level", "Total_SNPs",
        "Ancestral", "Derived", "Uncovered",
        "#ANC in par.", "#DER in par.", "#UNC in par.",
        "Total in par."
    ]

    # Check counts for haplogroup A
    rowA = df_out[df_out["Y-haplogroup"] == "A"].iloc[0]
    assert rowA["Total_SNPs"] == 2
    assert rowA["Ancestral"] == 1
    assert rowA["Derived"] == 1
    assert rowA["Uncovered"] == 0

    # Parental values come from mock
    assert rowA["Total in par."] == 30  # 3 × 10

def test_create_tree_basic():
    df = pd.DataFrame({
        "Y-haplogroup": ["A", "B"],
        "Level": [1, 0],
        "Derived": [2, 1],
        "Ancestral": [3, 4],
        "#DER in par.": [5, 0],
        "#ANC in par.": [6, 0],
        "Total_SNPs": [10, 20],
        "Total in par.": [11, 0]
    })

    chpar = {"A": "B", "B": None}

    m = mock_open()
    with patch("builtins.open", m):
        create_tree("A", level=1, df=df, chpar=chpar, file=m())

    written = "".join(call.args[0] for call in m().write.call_args_list)

    assert "|___> A" in written
    assert "|___> B" in written

def test_create_path_basic():
    chpar = {"A": "B", "B": "C", "C": None}

    path = create_path("A", df=None, chpar=chpar)

    assert path == ["C", "B", "A"]

def test_common_hap_basic():
    paths = [
        ["A", "B", "C", "D"],
        ["B", "C", "D"],
        ["C", "D"]
    ]

    assert common_hap(paths) == "D"

def test_recurse_leaf_node():
    node_info = {"children": {}, "samples": ["S1", "S2"]}
    avg_snps = {"A": 5}

    out = recurse("A", node_info, avg_snps)
    assert out == "'A [S1, S2]':5"

def test_recurse_internal_node():
    node_info = {
        "children": {
            "B": {"children": {}, "samples": []},
            "C": {"children": {}, "samples": ["S1"]}
        },
        "samples": []
    }
    avg_snps = {"A": 3, "B": 1, "C": 2}

    out = recurse("A", node_info, avg_snps)
    assert out == "('B':1,'C [S1]':2)A:3"

@patch("y_call.pulldown.recurse")
@patch("y_call.pulldown.pd.read_csv")
def test_generate_simple_grouped_tree_basic(mock_read_csv, mock_recurse):
    # Fake SNP lengths
    df_snps = pd.DataFrame({"node": ["A"], "length": [5]})
    mock_read_csv.return_value = df_snps

    # Fake recurse output
    mock_recurse.return_value = "(B:1)A:5"

    # Fake paths file
    fake_paths = "S1: A,B\n"

    m = mock_open()
    with patch("builtins.open", m):
        generate_simple_grouped_tree("paths.txt", "snps.txt", "out.nwk")

    # Check recurse was called with root node
    mock_recurse.assert_called_once()

    # Check Newick file written
    written = "".join(call.args[0] for call in m().write.call_args_list)
    assert written == "(B:1)A:5;"

@patch("y_call.pulldown.generate_simple_grouped_tree")
@patch("y_call.pulldown.Phylo.read")
@patch("y_call.pulldown.plt.figure")
def test_nwk_tree_basic(mock_fig, mock_read, mock_gen):
    mock_tree = MagicMock()
    mock_tree.get_terminals.return_value = []
    mock_tree.depths.return_value = {None: 1}
    mock_read.return_value = mock_tree

    mock_fig.return_value = MagicMock()

    nwk_tree("paths.txt", 10, 5, "snps.txt", "out.nwk", data=["S1", "S2"])

    mock_gen.assert_called_once()
    mock_read.assert_called_once()
    mock_fig.assert_called_once()

def test_unique_lineages_basic():
    df = pd.DataFrame({
        "Sample-ID": ["S1", "S2", "S3"],
        "Y-haplogroup": ["A", "B", "C"]
    })

    data = [
        "S1: A,B,C",
        "S2: B,C",
        "S3: C"
    ]

    # S1 shares A/B/C with others → not unique
    # S2 shares B/C → not unique
    # S3 shares C → not unique
    out = unique_lineages(df, data)
    assert out == ["S1"]

def test_unique_lineages_one_unique():
    df = pd.DataFrame({
        "Sample-ID": ["S1", "S2", "S3"],
        "Y-haplogroup": ["A", "B", "X"]
    })

    data = [
        "S1: A,B",
        "S2: B",
        "S3: X"
    ]

    # S3 has haplogroup X, not in any other path
    out = unique_lineages(df, data)
    assert out == ["S1", "S3"]

def test_get_mismatch_snps_basic():
    # Haplogroup tree: A -> B -> C (root)
    chpar = {"A": "B", "B": "C"}

    # Fake SNP table
    df_ch = pd.DataFrame({
        "Y-haplogroup": ["A", "A", "B", "C", "C"],
        "ref#": [5, 1, 3, 0, 2],
        "alt#": [0, 2, 1, 1, 0],
        "SNP-ID": ["rs1", "rs2", "rs3", "rs4", "rs5"]
    })

    # Expected mismatches:
    # A: rs1 (5 > 0)
    # B: rs3 (3 > 1)
    # C: rs5 (2 > 0)
    expected = ["rs1", "rs3", "rs5"]

    df_out = get_mismatch_snps("A", chpar, df_ch)

    assert sorted(df_out["SNP-ID"].tolist()) == expected

@patch("y_call.pulldown.nwk_tree")
@patch("y_call.pulldown.unique_lineages")
@patch("y_call.pulldown.create_path")
@patch("y_call.pulldown.pd.DataFrame.to_csv")
@patch("y_call.pulldown.pd.read_csv")
@patch("y_call.pulldown.pathlib.Path.mkdir")
@patch("y_call.pulldown.pathlib.Path.exists")
@patch("y_call.pulldown.os.path.exists")
def test_y_call_final_paths_and_summary(
    mock_os_exists, mock_path_exists, mock_mkdir,
    mock_read_csv, mock_to_csv, mock_create_path,
    mock_unique, mock_nwk
):
    """
    Tests the final block of y_call:
    - writing paths.txt
    - writing summary files
    - writing ages files
    - writing snps_branch.csv
    - computing unique lineages
    - optionally creating phylogeny
    """

    # -----------------------------
    # 1. Mock file existence logic
    # -----------------------------
    # os.path.exists("paths.txt") → True (file exists)
    mock_os_exists.return_value = True

    # Path.exists() calls:
    # - SNP file exists
    # - BED file exists
    # - BAM file exists
    # - ages file exists
    mock_path_exists.side_effect = [True, True, True, True]

    # -----------------------------
    # 2. Fake SNP dataset (OY)
    # -----------------------------
    df_snp = pd.DataFrame({
        "chrom": ["Y"],
        "pos": [100],
        "ref": ["A"],
        "alt": ["G"],
        "SNP-ID": ["rs1"],
        "Y-haplogroup": ["A"],
        "Level": [1]
    })

    # Fake BAM list
    df_bams = pd.DataFrame({
        "iid": ["S1"],
        "bam": ["file1.bam"],
        "sex": ["m"]
    })

    # Fake mismatch file
    df_mm = pd.DataFrame({"count": [0]})

    # Fake ages file
    df_ages = pd.DataFrame({
        "Y-haplogroup": ["A"],
        "Formed": [1000],
        "TMRCA": [800]
    })

    # Fake scores.csv (for unique_lineages)
    df_scores = pd.DataFrame({
        "Sample-ID": ["S1"],
        "Y-haplogroup": ["A"],
        "Level": [1],
        "#ANC in par./#DER in par.": [0.1],
        "Score": [5]
    })

    # pd.read_csv is called multiple times → return values in order
    mock_read_csv.side_effect = [
        df_snp,   # SNP dataset
        df_bams,  # BAM list
        df_mm,    # mismatch file
        df_ages,  # ages file
        df_scores # scores.csv for unique_lineages
    ]

    # -----------------------------
    # 3. Mock call_y_bam output
    # -----------------------------
    with patch("y_call.pulldown.call_y_bam") as mock_call_y:
        df_ch = pd.DataFrame({
            "Y-haplogroup": ["A"],
            "ref#": [5],
            "alt#": [0],
            "ref": ["A"],
            "alt": ["G"]
        })
        df_der = pd.DataFrame()

        mock_call_y.return_value = (df_ch, df_der)

        # -----------------------------
        # 4. Mock div_anc_der output
        # -----------------------------
        with patch("y_call.pulldown.div_anc_der") as mock_div:
            df_div = pd.DataFrame({
                "Y-haplogroup": ["A"],
                "Level": [1],
                "Total_SNPs": [10],
                "Ancestral": [3],
                "Derived": [7],
                "#ANC in par.": [1],
                "#DER in par.": [2]
            })
            mock_div.return_value = df_div

            # -----------------------------
            # 5. Mock create_path + unique_lineages
            # -----------------------------
            mock_create_path.return_value = ["C", "B", "A"]
            mock_unique.return_value = ["S1"]

            # -----------------------------
            # 6. Mock open()
            # -----------------------------
            m = mock_open(read_data="S1: C,B,A\n")
            with patch("builtins.open", m):

                # -----------------------------
                # 7. Run y_call
                # -----------------------------
                y_call(
                    bam_list="bam_list.tsv",
                    initial=0,
                    final=0,
                    base_qual=20,
                    map_qual=25,
                    database="testdb",
                    reference_genome="hg19",
                    create_phylogeny="Y",
                    width=10,
                    height=10,
                    transitions="N",
                    translation="translation.csv",
                    ex_limit=5,
                    ages="ages.csv"
                )

    # -----------------------------
    # 8. Assertions
    # -----------------------------

    # Path writing must have occurred
    assert m().write.called

    # Summary file must have been written
    assert mock_to_csv.call_count >= 3

    # Unique lineages must have been computed
    mock_unique.assert_called_once()

    # Phylogeny must have been created
    mock_nwk.assert_called_once()

    # create_path must have been used
    mock_create_path.assert_called()
