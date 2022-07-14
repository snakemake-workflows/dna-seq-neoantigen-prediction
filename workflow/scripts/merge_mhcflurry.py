import sys

sys.stderr = open(snakemake.log[0], "w")

import os
import pandas as pd
import numpy as np


## highlight the difference between mutated neopeptide and wildtype
def diffEpitope(e1,e2):
    if str(e2) == 'nan':
        return(e1)
    e1 = str(e1)
    e2 = str(e2)
    diff_pos = [i for i in range(len(e1)) if e1[i] != e2[i]]
    e_new = e1
    e2_new = e2
    for p in diff_pos:
        e_new = e_new[:p] + e_new[p].lower() + e_new[p+1:]
        e2_new = e2_new[:p] + e2_new[p].lower() + e2_new[p+1:]
    return(e_new)

info = pd.read_csv(snakemake.input[0])
tumor = pd.read_csv(snakemake.input[1])
normal = pd.read_csv(snakemake.input[2])
outfile = snakemake.output[0]


tumor = tumor[["source_sequence_name","peptide","allele","affinity","percentile_rank"]]
tumor = tumor.pivot_table(["affinity","percentile_rank"],["source_sequence_name","peptide"],"allele").reset_index()
tumor.columns = tumor.columns.map("-".join)
tumor = tumor.rename(columns={col: col.replace("-","") for col in tumor.columns if col.endswith("-")})

normal = normal[["source_sequence_name","peptide","allele","affinity","percentile_rank"]]
normal = normal.pivot_table(["affinity","percentile_rank"],["source_sequence_name","peptide"],"allele").reset_index()
normal.columns = normal.columns.map("-".join)
normal = normal.rename(columns={col: col.replace("-","") for col in normal.columns if col.endswith("-")})

merged = tumor.merge(normal, how="left", on=["source_sequence_name"])

merged = merged.rename(columns={col: col.replace("_y","_normal") for col in merged.columns}).rename(columns={col: col.replace("_x","_tumor") for col in merged.columns})
## add info
info = info.rename(columns={"id":"ID","gene_id":"Gene_ID","gene_name":"Gene_Symbol","strand":"Strand","positions":"Variant_Position","chrom":"Chromosome","somatic_aa_change":"Somatic_AminoAcid_Change"})
merged_dataframe = merged.merge(info, how="left", left_on="source_sequence_name", right_on="ID")

merged_dataframe["peptide_tumor"]=merged_dataframe[["peptide_tumor","peptide_normal"]].apply(lambda x: diffEpitope(*x), axis=1)
## Are all possible variants in the peptide ("Cis") or not ("Trans")
merged_dataframe["Variant_Orientation"] = "Cis"
trans = merged_dataframe.nvariant_sites > merged_dataframe.nvar
merged_dataframe.loc[trans, "Variant_Orientation"] = "Trans"

## check misssense/silent mutation status
nonsilent = merged_dataframe.peptide_tumor != merged_dataframe.peptide_normal
merged_dataframe = merged_dataframe[nonsilent]
data = merged_dataframe.drop_duplicates(subset=["Gene_ID","offset","peptide_tumor","Somatic_AminoAcid_Change"])

### Delete Stop-Codon including peptides
data = data[data.peptide_tumor.str.count("x") == 0]
data = data[data.peptide_tumor.str.count("X") == 0]
data.to_csv(outfile, index=False, sep = '\t')





