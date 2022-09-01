import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl

columns_mapping = {
    "gene_name": "gene",
    "normal_peptide": "mutation.wildTypeXmer",
    "tumor_peptide": "mutation.mutatedXmer",
    "freq": "dnaVariantAlleleFrequency",
}

candidates = (
    pl.read_tsv(snakemake.input.microphaser, sep="\t", quote="")
    .rename(columns_mapping)
    .with_column(pl.lit(snakemake.wildcards.group).alias("patientIdentifier"))
)

candidates.write_csv(snakemake.output.neo_fox, sep="\t", quote="")
