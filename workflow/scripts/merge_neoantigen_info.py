import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

def get_minimum_rank_per_peptide(df: pd.DataFrame, rank_type: str) -> pd.DataFrame:
    df = df.set_index(['id', 'pos_in_id_seq', 'pep_seq', 'allele'])
    rank_col = f"{rank_type}_rank"
    score_col = f"{rank_type}_score"
    prefix = f"top_{rank_col}_"
    columns_to_keep = ['bind_core', rank_col, score_col]
    return df.loc[ df.groupby(['pep_seq', 'id'])[rank_col].idxmin(), columns_to_keep ].reset_index(level='allele').sort_index().drop_duplicates().add_prefix(prefix)

def select_columns(mhc: pd.DataFrame) -> pd.DataFrame:
    rank_cols = [c for c in mhc.columns if "Rank" in c]
    affinity_cols = [c for c in mhc.columns if "nM" in c]
    mhc_cols = ["Pos", "ID", "Peptide"] + rank_cols + affinity_cols + ["NB"]
    mhc = mhc[mhc_cols]
    mhc["Rank_min"] = mhc[rank_cols].min(axis=1)
    mhc["Aff_min"] = mhc[affinity_cols].min(axis=1)
    mhc["Top_rank_HLA"] = mhc[rank_cols].idxmin(axis=1)
    mhc["Top_affinity_HLA"] = mhc[affinity_cols].idxmin(axis=1)
    mhc["Top_rank_HLA"] = mhc["Top_rank_HLA"].str.replace("_Rank","")
    mhc["Top_affinity_HLA"] = mhc["Top_affinity_HLA"].str.replace("_nM","")
    return mhc

def get_filtered_per_sample(sample: pd.DataFrame, alias: str) -> pd.DataFrame:
    common_info = sample.set_index(['id', 'pos_in_id_seq', 'pep_seq']).loc[:, ['ave_el_score', 'num_binders']].reset_index(level=['id', 'pos_in_id_seq']).drop_duplicates().set_index(['id', 'pos_in_id_seq'], append=True)
    sample_el = get_minimum_rank_per_peptide(sample, "el")
    sample_ba = get_minimum_rank_per_peptide(sample, "ba")
    sample_filtered = sample_el.join(sample_ba)
    return sample_filtered.join(common_info, how='left').assign(alias=alias).set_index('alias', append=True)

def tidy_info(info: pd.DataFrame, alias: str) -> pd.DataFrame:
    info = info.set_index(['id', 'transcript', 'gene_id', 'gene_name', 'chrom', 'offset', 'frame', 'freq', 'credible_interval', 'depth', 'strand'])
    int_cols = ['nvar', 'nsomatic', 'nvariant_sites', 'nsomvariant_sites']
    info[int_cols] = info[int_cols].astype('int32')
    num_var_tidy = info.assign(ngermline = lambda x: x.nvar - x.nsomatic).melt(var_name='alias', value_name='num_var', value_vars=['ngermline', 'nsomatic'], ignore_index=False).replace({'ngermline': 'normal', 'nsomatic': alias}).set_index('alias', append=True)
    num_var_sites_tidy = info.assign(ngermvariant_sites = lambda x: x.nvariant_sites - x.nsomvariant_sites).melt(var_name='alias', value_name='num_var_sites', value_vars=['ngermvariant_sites', 'nsomvariant_sites'], ignore_index=False).replace({'ngermvariant_sites': 'normal', 'nsomvariant_sites': alias}).set_index('alias', append=True)
    sites_tidy = info.copy()
    sites_tidy['sites'] = sites_tidy['variant_sites'].apply(lambda x: set(str(x).split('|')))
    sites_tidy['somatic_sites'] = sites_tidy['somatic_positions'].apply(lambda x: set(str(x).split('|')))
    sites_tidy['germline_sites'] = sites_tidy[['sites', 'somatic_sites']].apply(lambda row: [s for s in row['sites'] if s not in row['somatic_sites']], axis=1)
    sites_tidy = sites_tidy.melt(var_name='alias', value_name='genomic_pos', value_vars=['somatic_sites', 'germline_sites'], ignore_index=False).replace({'somatic_sites': 'normal', 'germline_sites': 'tumor_alias'}).set_index('alias', append=True)['genomic_pos'].apply(lambda r: '|'.join(r))
    sites_tidy = info.melt(var_name='alias', value_name='genomic_pos', value_vars=['somatic_positions', 'germline_positions'], ignore_index=False).replace({'somatic_positions': 'tumor_resection', 'germline_positions': 'normal'}).set_index('alias', append=True)

    seq_tidy = info.melt(var_name='alias', value_name='sequence', value_vars=['normal_sequence', 'mutant_sequence'], ignore_index=False).replace({'normal_sequence': 'normal', 'mutant_sequence': alias}).set_index('alias', append=True)

def merge_data_frames(info: pd.DataFrame, tumor: pd.DataFrame, normal: pd.DataFrame) -> pd.DataFrame:
    # get and merge tumor and normal
    tumor_filtered = get_filtered_per_sample(tumor, snakemake.wildcards.tumor_alias)
    normal_filtered = get_filtered_per_sample(normal, "normal")
    all_filtered = pd.concat([tumor_filtered, normal_filtered]).reset_index(level='pep_seq').sort_index()
    info_tidy = tidy_info(info, snakemake.wildcards.tumor_alias)

    # tidy up info
    all_annotated = all_filtered.merge(info, on='id', how='left')
#    tumor = select_columns(tumor)
#    normal = select_columns(normal)
#    id_length = len(tumor.ID[0])
#    print(info.columns)
#    info["ID"] = info["id"].astype(str).str[:id_length]
#    merged_mhc = tumor.merge(normal,how='left', on=['Pos','ID'])
#    merged_mhc = merged_mhc
#        .rename(
#            columns={col: col.replace("_y","_normal") for col in merged_mhc.columns}
#            )
#        .rename(
#            columns={col: col.replace("_x","_tumor") for col in merged_mhc.columns}
#            )
#    info = info.rename(
#        columns={
#            "gene_id":"Gene_ID",
#            "gene_name":"Gene_Symbol",
#            "strand":"Strand",
#            "positions":"Variant_Position",
#            "chrom":"Chromosome",
#            "somatic_aa_change":
#            "Somatic_AminoAcid_Change"
#        })
#    merged_dataframe = merged_mhc.merge(info, how='left', on = 'ID')

    merged_dataframe["Peptide_tumor"] = merged_dataframe[["Peptide_tumor","Peptide_normal"]].apply(lambda x: diffEpitope(*x), axis=1)
    ## Are all possible variants in the peptide ("Cis") or not ("Trans")
    merged_dataframe["Variant_Orientation"] = "Cis"
    trans = merged_dataframe.nvariant_sites > merged_dataframe.nvar
    merged_dataframe.loc[trans, "Variant_Orientation"] = "Trans"

    ## check misssense/silent mutation status
    nonsilent = merged_dataframe.Peptide_tumor != merged_dataframe.Peptide_normal
    merged_dataframe = merged_dataframe[nonsilent]
    merged_dataframe = merged_dataframe.drop_duplicates(subset=["transcript","offset","Peptide_tumor","Somatic_AminoAcid_Change"])

    data = merged_dataframe[["ID","transcript","Gene_ID","Gene_Symbol","Chromosome","offset","freq","depth",
    "Somatic_AminoAcid_Change", "nvar", "nsomatic", "somatic_positions", "Peptide_tumor","NB_tumor","Rank_min_tumor","Aff_min_tumor",
    "Top_rank_HLA_tumor","Top_affinity_HLA_tumor","Peptide_normal","NB_normal",
    "Rank_min_normal","Aff_min_normal","Top_rank_HLA_normal","Top_affinity_HLA_normal"]]

    data.columns = ["ID","Transcript_ID","Gene_ID","Gene_Symbol","Chromosome","Position","Frequency","Read_Depth",
    "Somatic_AminoAcid_Change", "nvar", "nsomatic", "somatic_positions", "Peptide_tumor","BindingHLAs_tumor","Rank_min_tumor","Affinity_min_tumor",
    "Top_rank_HLA_tumor","Top_affinity_HLA_tumor","Peptide_normal","BindingHLAs_normal",
    "Rank_min_normal","Aff_min_normal","Top_rank_HLA_normal","Top_affinity_HLA_normal"]

    # data = data[data.BindingHLAs_tumor > 0]
    # data = data[(data.NB_normal.isna()) | (data.NB_normal == 0)]
    #data = data[(data.BindingHLAs_normal == 0)]

    ### Delete Stop-Codon including peptides
    data = data[data.Peptide_tumor.str.count("x") == 0]
    data = data[data.Peptide_tumor.str.count("X") == 0]
    data.sort_values(["Chromosome", "somatic_positions"], inplace=True)
    ### Remove Duplicate kmers
    data = data.drop_duplicates(["Transcript_ID", "Peptide_tumor", "Somatic_AminoAcid_Change", "Peptide_normal"])

    return data


## highlight the difference between mutated neopeptide and wildtype
def diffEpitope(e1,e2):
    if str(e2) == 'nan' or str(e2) == '':
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
    info = pd.read_csv(snakemake.input.info, sep = '\t', dtype=str)
    tumor = pd.read_csv(snakemake.input.neo, sep = '\t', dtype={'pos_in_id_seq': str})
    normal = pd.read_csv(snakemake.input.normal, sep = '\t', dtype={'pos_in_id_seq': str})
    data = merge_data_frames(info, tumor, normal)
    data.to_csv(snakemake.output[0], index=False, sep = '\t')

if __name__ == '__main__':
    sys.exit(main())
