import pandas as pd

## load data table
data = pd.read_csv(snakemake.input["table"], sep='\t')

## Merge transcript count
transcript_count = pd.read_csv(snakemake.input["counts"], sep='\t')
transcript_count = transcript_count[["target_id", "tpm"]]
transcript_count.columns = ["Transcript_ID", "TPM"]
data = data.merge(transcript_count, on="Transcript_ID", how="left")

data.to_csv(snakemake.output[0], sep='\t', index=False)
