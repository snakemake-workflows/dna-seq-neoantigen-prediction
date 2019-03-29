import sys
import os
import pandas as pd
import numpy as np


def merge(info, tumor, normal, outfile):
    rank_cols = [c for c in tumor.columns if "Rank" in c]
    affinity_cols = [c for c in tumor.columns if "nM" in c]
    mhc_cols = ["ID"] + ["Peptide"] + rank_cols + affinity_cols + ["NB"]
    tumor = tumor[mhc_cols]
	normal = normal[mhc_cols]
	for mhc in [tumor, normal]
		mhc["Rank_min"] = mhc[rank_cols].min(axis=1)
		mhc["Aff_min"] = mhc[affinity_cols].min(axis=1)
		mhc["Top_rank_HLA"] = mhc[rank_cols].idxmin(axis=1)
		mhc["Top_affinity_HLA"] = mhc[affinity_cols].idxmin(axis=1)
		mhc["Top_rank_HLA"] = mhc["Top_rank_HLA"].str.replace("Rank_","")
		mhc["Top_affinity_HLA"] = mhc["Top_affinity_HLA"].str.replace("nM_","")
    info["ID"] = info["id"].astype(str).str[:-1]

    merged_mhc = tumor.merge(normal,how='left', on='ID')
	merged_mhc = merged_mhc.rename(columns={col: col.replace("_y","_normal") for col in merged_mhc.columns}).rename(columns={col: col.replace("_x","_tumor") for col in merged_mhc.columns})

	merge_mhc = merge_mhc.rename(columns={"gene_id":"Gene_ID","gene_name":"Gene_Symbol","strand":"Strand","positions":"Variant_Position","chrom":"Chromosome","aa_change":"AminoAcid_Change"})

    
	merged_dataframe = merged_mhc.merge(info, how='left', on = 'ID')



    merge_mhc["Peptide_tumor"]=merge_mhc[["Peptide_tumor","Peptide_normal"]].apply(lambda x: diffEpitope(*x), axis=1)
    ## Are all possible variants in the peptide ("Cis") or not ("Trans")
    merge_mhc["Variant_Orientation"] = "Cis"
    trans = merge_mhc.nvariant_sites_tumor > merge_mhc.nvar_tumor
    merge_mhc.loc[trans, "Variant_Orientation"] = "Trans"

    ## check misssense/silent mutation status
    nonsilent = merge_mhc.Peptide_tumor != merge_mhc.Peptide_normal
    merge_mhc = merge_mhc[nonsilent]
    merge_mhc = merge_mhc.drop_duplicates(subset=["Gene_ID","offset","Peptide_tumor","AminoAcid_Change"])

    merge_mhc.to_csv(outfile, index=False, sep = '\t')

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


def main():
    info = pd.read_csv(snakemake.input[0], sep = '\t')
    tumor = pd.read_csv(snakemake.input[1], sep = '\t')
	normal = pd.read_csv(snakemake.input[2], sep = '\t')
    outfile = snakemake.output[0]
    merge(info, tumor, normal, outfile)

if __name__ == '__main__':
    sys.exit(main())
