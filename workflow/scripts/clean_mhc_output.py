import sys

sys.stderr = open(snakemake.log[0], "w")

import os
import pandas as pd

# assumptions of this script about netMHCpan or netMHCIIpan:
# * version 4.1
# * output generated via `-xls` option
# * generated with the `-BA` option to include binding affinity prediction

# The mapping of index column names used here to original names in netMHCpan files is:
INDEX_NAMES = {
    "Pos": "position_in_protein_sequence",
    "Peptide": "peptide_sequence",
    "ID": "peptide_ID",
    "Ave": "average_el_score",
    "NB": "number_of_binders",
}
# The mapping of column names used here to original names in netMHCpan files is:
COLUMN_NAMES = {
    "BA-score": "binding_affinity_score",
    "BA_Rank": "binding_affinity_percent_rank",
    "EL-score": "elution_ligang_score",
    "EL_Rank": "elution_ligand_percent_rank",
    "core": "binding_core",
    "icore": "interaction_core",
}


def parse_file(mhc_in: FileIO):
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

    # It's a compound header over two rows and a compound row index in the initial
    # three and final two columns of the table. For some reason, the final two
    # columns are added to the index but not removed from the table, so we do this
    # manually with `.iloc[]``.
    data = pd.read_csv(
        mhc_in, sep="\t", header=[0, 1], index_col=[0, 1, 2, -2, -1]
    ).iloc[:, :-2]

    # With two lines of header parsed into a MultiIndex, pandas only uses the
    # first column name in index_col as an entry. Obviously the following code
    # assumes that these are the first three and last two columns of the data file.
    data.index.names = list(INDEX_NAMES.values())

    # Entries of the header MultiIndex need to be fixed, there doesn't seem to be
    # any way to automatically do this during read_csv.
    cols = pd.DataFrame(data.columns.to_list(), columns=["allele", "info"])

    # fix up the columns and reassign
    cols.loc[cols["allele"].str.endswith("_level_0"), "allele"] = pd.NA
    cols = cols.fillna(method="ffill")
    data.columns = pd.MultiIndex.from_frame(cols)

    # Turn into longer table with one HLA Allele per row instead of MultiIndex
    # header, rename columns to something readable and turn index into columns.
    data = data.stack(level="allele").rename(columns=COLUMN_NAMES).reset_index()

    return data


all_data = pd.concat((parse_file(f) for f in snakemake.input), axis="index")

all_data.to_csv(snakemake.output.joined_mhc_out, sep="\t", index=False)
