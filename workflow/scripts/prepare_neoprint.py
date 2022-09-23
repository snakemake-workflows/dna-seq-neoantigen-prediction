import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from typing import Tuple

def highlight_peptides_diff(tumor_p: str, normal_p: str) -> Tuple[str, str]:
    """
    Highlight the difference between mutated neopeptide and normal peptide
    """
    if normal_p == "nan" or normal_p == "NA" or normal_p == "":
        return (tumor_p, normal_p)
    assert len(tumor_p) == len(
        normal_p
    ), f"Tumor peptide '{tumor_p}' and normal peptide '{normal_p}' have different lengths."
    diff_pos = [i for i in range(len(tumor_p)) if tumor_p[i] != normal_p[i]]
    tp_changed = tumor_p
    np_changed = normal_p
    for p in diff_pos:
        tp_changed = tp_changed[:p] + tp_changed[p].lower() + tp_changed[p + 1 :]
        np_changed = np_changed[:p] + np_changed[p].lower() + np_changed[p + 1 :]
    return (tp_changed, np_changed)

all_neopeptides = pd.read_csv(snakemake.input.neopeptides, sep="\t")

# If we leave in any dots, datavzrd will interpret this as attributes
all_neopeptides.columns = all_neopeptides.columns.str.replace(".", "_", regex=False)

# Aggregate multiple identical entries that differ only in 'id' and 'transcript'
# into one, taking the first 'id' and collecting all 'transcript's into a '|'-separated
# list.
# TODO: Remove this redundancy from microphaser output before passing it along to other
# tools.
cols = [c for c in all_neopeptides.columns if c not in ["id", "transcript"]]
aggregation_functions = {
    "id": lambda i: list(i)[0],
    "transcript": lambda t: "|".join(list(t)),
}
all_neopeptides = (
    all_neopeptides.groupby(cols, dropna=False)
    .agg(aggregation_functions)
    .reset_index()
    .explode("id")
)

# highlight mutations in the original Xmers
all_neopeptides[["mutation_mutatedXmer", "mutation_wildTypeXmer"]] = (
    pd.DataFrame(
        all_neopeptides
        .fillna({"mutation_wildTypeXmer": ""})
        .apply(lambda row: highlight_peptides_diff(row["mutation_mutatedXmer"], row["mutation_wildTypeXmer"]), axis="columns")
        .tolist()
    )
)

# create new purity-adjusted DNA variant allele frequency
all_neopeptides["purity_adjusted_DNA_VAF"] = all_neopeptides["dnaVariantAlleleFrequency"] / snakemake.params.purity

# round all floats to the specified decimals
all_neopeptides = all_neopeptides.round(decimals=5)

# define important columns to move to the left of the table

important_cols_general = [
            "gene",
            "mutation_mutatedXmer",
            "mutation_wildTypeXmer",
            "purity_adjusted_DNA_VAF",
            "imputedGeneExpression",
        ]

important_cols_one = [
            "PRIME_best_rank",
            "PRIME_best_score",
            "PRIME_best_peptide",
            "PRIME_best_allele",
            "Best_rank_MHCI_9mer_score",
            "Best_rank_MHCI_9mer_score_WT",
            "Best_rank_MHCI_9mer_epitope",
            "Best_rank_MHCI_9mer_epitope_WT",
            "Best_rank_MHCI_9mer_allele",
            "Best_rank_MHCI_9mer_allele_WT",
        ]

important_cols_two = [
            "MixMHC2pred_best_rank",
            "MixMHC2pred_best_peptide",
            "MixMHC2pred_best_allele",
            "Best_rank_MHCII_score",
            "Best_rank_MHCII_score_WT",
            "Best_rank_MHCII_score_epitope",
            "Best_rank_MHCII_score_epitope_WT",
            "Best_rank_MHCII_score_allele",
            "Best_rank_MHCII_score_allele_WT",
        ]

important_cols = important_cols_general + important_cols_one + important_cols_two


mhc_one = (
    all_neopeptides[ important_cols + [ col for col in all_neopeptides.columns if col not in important_cols ] ]
    .sort_values(by = ["PRIME_best_rank", "MixMHC2pred_best_rank"])
    .groupby("PRIME_best_rank")
    .head(n=1)
)

mhc_one.to_csv(snakemake.output.mhc_one, sep="\t", index=False)

# move important columns to the left of the table
mhc_two = (
    all_neopeptides[ important_cols_general + important_cols_two + important_cols_one + [ col for col in all_neopeptides.columns if col not in important_cols ] ]
    .sort_values(by = ["MixMHC2pred_best_rank", "PRIME_best_rank"])
    .groupby("MixMHC2pred_best_rank")
    .head(n=1)
)

mhc_two.to_csv(snakemake.output.mhc_two, sep="\t", index=False)
