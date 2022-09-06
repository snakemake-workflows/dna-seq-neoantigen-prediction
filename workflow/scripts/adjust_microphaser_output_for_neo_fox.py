import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

columns_mapping = {
    "gene_name": "gene",
    "normal_peptide": "mutation.wildTypeXmer",
    "tumor_peptide": "mutation.mutatedXmer",
    "freq": "dnaVariantAlleleFrequency",
}

candidates = (
    pd.read_csv(snakemake.input.microphaser, sep="\t", quoting=3)
    .rename(columns=columns_mapping)
    .assign(patientIdentifier=snakemake.wildcards.group)
)

candidates.to_csv(snakemake.output.neo_fox, sep="\t", quoting=3)
