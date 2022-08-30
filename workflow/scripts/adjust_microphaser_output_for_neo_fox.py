import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl

candidates = pl.read_tsv(snakemake.input.candidates, sep="\t")