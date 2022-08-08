import sys

sys.stderr = open(snakemake.log[0], "w")

import os
import pandas as pd

from itertools import cycle

# assumptions of this script about netMHCpan or netMHCIIpan:
# * version 4.1
# * output generated via `-xls` option
# * generated with the `-BA` option to include binding affinity prediction

# The mapping of index column names used here to original names in netMHCpan files
# is (please excuse the pd.NA tuples, they make header and index handling 
# easier further down the line):
INDEX_NAMES = {
    (pd.NA, "Pos"): "position_in_protein_sequence",
    (pd.NA, "Peptide"): "peptide_sequence",
    (pd.NA, "ID"): "peptide_ID",
    (pd.NA, "Ave"): "average_el_score",
    (pd.NA, "NB"): "number_of_binders",
}

if snakemake.wildcards.mhc == "net_mhc_pan":
    # The mapping of column names used here to original names in netMHCpan files is:
    COLUMN_NAMES = {
        "BA-score": "binding_affinity_score",
        "BA_Rank": "binding_affinity_percent_rank",
        "EL-score": "elution_ligang_score",
        "EL_Rank": "elution_ligand_percent_rank",
        "core": "binding_core",
    }
elif snakemake.wildcards.mhc == "net_mhc_two_pan":
    # The mapping of column names used here to original names in netMHCIIpan files is:
    COLUMN_NAMES = {
        "Score_BA": "binding_affinity_score",
        "Rank_BA": "binding_affinity_percent_rank",
        "Score": "elution_ligang_score",
        "Rank": "elution_ligand_percent_rank",
        "Core": "binding_core",
    }
else:
    sys.exit(f"Wildcard `mhc` has unknown value: {snakemake.wildcards.mhc}")

COLUMNS_TO_DROP = {
    # I have not found any docs or indication in the manuscript, what the column
    # `Target` is. It might be the column `Exp_Bind` from Figure 1B here (it's
    # the only column that only appears in netMHCIIpan output and not netMHCpan
    # output):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7319546/figure/F1/
    # If so, it is only for benchmarking purposes according to the docs. And it
    # is always NA, even in their manuscript. So we simply remove it here.
    "Target",
    # There doesn't seem to be an `icore` equivalent in netMHCIIpan output.
    "icore",
    # There doesn't seem to be an `nM` equivalent in netMHCpan output.
    "nM",
}


def parse_file(mhc_in: str):
    """
    Parse an netMHCpan or netMHCIIpan output file from the `-xls -xlsfile <filename>`
    directive into a tidy pandas data frame.
    """
    if os.path.getsize(mhc_in) == 0:
        # Short-circuit empty files, but generate correct header.
        return pd.DataFrame(
            columns=list(COLUMN_NAMES.values())
            + ["allele"]
            + list(INDEX_NAMES.values())
        )

    # We need to fix the utterly broken headers first

    # parse first header line into pandas.Series and name it
    first_header_line = pd.read_csv(mhc_in, nrows=1, header=None, sep="\t").iloc[0, :]
    first_header_line.name = "allele"

    # parse second header line into pandas.Series and name it
    second_header_line = pd.read_csv(
        mhc_in, skiprows=1, nrows=1, header=None, sep="\t"
    ).iloc[0, :]
    second_header_line.name = "column_name"

    header = pd.concat([first_header_line, second_header_line], axis="columns")
    header = header.fillna(method="ffill")
    header.loc[
        header.column_name.isin({"Pos", "Peptide", "ID", "Target", "Ave", "NB"}),
        "allele",
    ] = pd.NA

    # It's a compound header over two rows and a compound row index in the initial
    # three and final two columns of the table. For some reason, the final two
    # columns are added to the index but not removed from the table, so we do this
    # manually with `.iloc[]``.
    data = pd.read_csv(mhc_in, sep="\t", skiprows=2, header=None)

    data.columns = pd.MultiIndex.from_frame(header)

    # remove columns only present in one of the two tools' output
    columns_to_keep = [
        col
        for col in list(data.columns.get_level_values("column_name"))
        if col not in COLUMNS_TO_DROP
    ]
    idx = pd.IndexSlice
    data = data.loc[:, idx[:, columns_to_keep]]

    # properly set row index columns (values to repeat in every row while doing the stack below)
    data = data.set_index(list(INDEX_NAMES.keys()))
    data.index.set_names(INDEX_NAMES, inplace=True)

    # Turn into longer table with one HLA Allele per row instead of MultiIndex
    # header, rename columns to something readable and turn index into columns.
    data = data.stack(level="allele").rename(columns=COLUMN_NAMES).reset_index()

    return data


all_data = pd.concat((parse_file(f) for f in snakemake.input), axis="index")

all_data.to_csv(snakemake.output.joined_mhc_out, sep="\t", index=False)
