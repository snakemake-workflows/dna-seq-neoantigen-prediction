import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

data = pd.read_csv(snakemake.input.tsv, sep="\t")
data.to_excel(snakemake.output.xlsx, index=False)