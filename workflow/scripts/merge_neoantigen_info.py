import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from typing import List, Tuple


def get_best_rank_per_peptide(df: pd.DataFrame, rank_type: str) -> pd.DataFrame:
    df = df.set_index(["id", "pos_in_id_seq", "pep_seq", "allele"])
    rank_col = f"{rank_type}_rank"
    score_col = f"{rank_type}_score"
    prefix = f"top_{rank_col}_"
    columns_to_keep = ["bind_core", rank_col, score_col]
    return (
        df.loc[df.groupby(["pep_seq", "id"])[rank_col].idxmin(), columns_to_keep]
        .reset_index(level="allele")
        .add_prefix(prefix)
        .reset_index()
        .drop_duplicates()
    )


def get_filtered_per_alias(sample: pd.DataFrame, alias: str) -> pd.DataFrame:
    merge_cols = ["id", "pos_in_id_seq", "pep_seq"]
    common_info = (
        sample.loc[:, merge_cols + ["ave_el_score", "num_binders"]]
        .drop_duplicates()
    )
    sample_el = get_best_rank_per_peptide(sample, "el")
    sample_ba = get_best_rank_per_peptide(sample, "ba")
    sample_filtered = sample_el.merge(sample_ba, on=merge_cols)
    return (
        sample_filtered.merge(common_info, how="left", on=merge_cols)
        .assign(alias=alias)
    )


def highlight_peptides_diff(tumor_p: str, normal_p: str) -> Tuple[str, str]:
    """
    Highlight the difference between mutated neopeptide and normal peptide
    """
    if normal_p == "nan" or normal_p == "":
        return (tumor_p, normal_p)
    assert len(tumor_p) == len(normal_p), f"Tumor peptide '{tumor_p}' and normal peptide '{normal_p}' have different lengths."
    diff_pos = [i for i in range(len(tumor_p)) if tumor_p[i] != normal_p[i]]
    tp_changed = tumor_p
    np_changed = normal_p
    for p in diff_pos:
        tp_changed = tp_changed[:p] + tp_changed[p].lower() + tp_changed[p + 1 :]
        np_changed = np_changed[:p] + np_changed[p].lower() + np_changed[p + 1 :]
    return (tp_changed, np_changed)


def diff_tumor_normal_peptides(
    group: pd.DataFrame, column: str, tumor_alias: str
) -> pd.DataFrame:
    group = group
    normal_pep = group.loc[group["alias"] == "normal", column].fillna("")
    if normal_pep.empty:
        normal_pep = ""
    else:
        normal_pep = normal_pep.squeeze()
    tumor_pep = group.loc[group["alias"] == tumor_alias, column].fillna("").squeeze()
    # Silent mutations should not be included in microphaser output.
    if normal_pep == tumor_pep:
        raise ValueError(
            f"For peptide '{group['id'][0]}' the normal and the tumor peptide have an identical sequence ({normal_pep}).\n"
            "Please fix this upstream or comment out this check to ignore this problem.\n"
        )
    # Remove groups where the tumor peptide contains a stop codon.
    # TODO: Maybe this should be a hard fail complaining to fix this upstream?
    # TODO: Write out warning.
    if "X" in tumor_pep:
        print(f"Warning: ", file=sys.stderr)
        return group.loc[[], :]
    t_diff, n_diff = highlight_peptides_diff(tumor_pep, normal_pep)
    group.loc[group["alias"] == tumor_alias, column] = t_diff
    group.loc[group["alias"] == "normal", column] = n_diff
    return group


def tidy_info(info: pd.DataFrame, tumor_alias: str) -> pd.DataFrame:
    """
    Get the -o info output of the microphaser filter command into tidy data format.
    """
    info = info.rename(
        columns={"credible_interval": "freq_credible_interval"}
    ).set_index(
        [
            "id",
            "transcript",
            "gene_id",
            "gene_name",
            "chrom",
            "offset",
            "frame",
            "freq",
            "freq_credible_interval",
            "depth",
            "strand",
        ]
    )
    int_cols = ["nvar", "nsomatic", "nvariant_sites", "nsomvariant_sites"]
    info[int_cols] = info[int_cols].astype("int32")
    # TODO: Ensure that microphaser output contains only one entry per id.
    # If there is more than one entry per index, ensure that they are identical
    if len(info.groupby(info.index).filter(lambda g: (g.nunique() > 1).any())) > 0:
        raise ValueError(
            f"Found multiple differing entries for an 'id' in file '{snakemake.input.info}'. Please ensure that entries are unique per 'id'.\n"
        )
    # Always take the first entry for each index.
    info = info.groupby(info.index).head(1)
    # TODO: Ensure that microphaser output is tidy data, with one row each for tumor and normal.
    # TODO: factor out tidying of column pairs into a function.
    num_var_in_pep_tidy = (
        info.assign(ngermline=lambda x: x.nvar - x.nsomatic)
        .melt(
            var_name="alias",
            value_name="num_var_in_pep",
            value_vars=["ngermline", "nsomatic"],
            ignore_index=False,
        )
        .replace({"ngermline": "normal", "nsomatic": tumor_alias})
        .set_index("alias", append=True)
    )
    num_var_sites_tidy = (
        info.assign(ngermvariant_sites=lambda x: x.nvariant_sites - x.nsomvariant_sites)
        .melt(
            var_name="alias",
            value_name="num_var_sites",
            value_vars=["ngermvariant_sites", "nsomvariant_sites"],
            ignore_index=False,
        )
        .replace({"ngermvariant_sites": "normal", "nsomvariant_sites": tumor_alias})
        .set_index("alias", append=True)
    )
    genomic_pos_tidy = (
        info.melt(
            var_name="alias",
            value_name="genomic_pos",
            value_vars=["somatic_positions", "germline_positions"],
            ignore_index=False,
        )
        .replace({"somatic_positions": tumor_alias, "germline_positions": "normal"})
        .set_index("alias", append=True)
    )
    aa_changes_tidy = (
        info.melt(
            var_name="alias",
            value_name="aa_changes",
            value_vars=["somatic_aa_change", "germline_aa_change"],
            ignore_index=False,
        )
        .replace({"somatic_aa_change": tumor_alias, "germline_aa_change": "normal"})
        .set_index("alias", append=True)
    )
    nt_seq_tidy = (
        info.melt(
            var_name="alias",
            value_name="nt_seq",
            value_vars=["normal_sequence", "mutant_sequence"],
            ignore_index=False,
        )
        .replace({"normal_sequence": "normal", "mutant_sequence": tumor_alias})
        .set_index("alias", append=True)
    )
    all_tidy = num_var_in_pep_tidy.join(
        [num_var_sites_tidy, genomic_pos_tidy, aa_changes_tidy, nt_seq_tidy]
    )
    return all_tidy.reset_index()


def check_duplicates(df: pd.DataFrame, cols: List[str], specific_error: str = ""):
    if (
        sum(
            df.duplicated(
                subset=cols
            )
        )
        > 0
    ):
        duplicates = all_annotated[
            df.duplicated(
                subset=cols
            )
        ]
        cols_str = '", "'.join(cols)
        raise ValueError(
            f'Found multiple rows with identical [ "{cols_str}" ] entries.\n'
            "This indicates an upstream issue, please fix this.\n"
            f"{specific_error}"
            "The offending entries are:\n"
            f"{duplicates}\n"
        )


def merge_data_frames(
    info: pd.DataFrame, tumor: pd.DataFrame, normal: pd.DataFrame, tumor_alias: str
) -> pd.DataFrame:
    # get and merge tumor and normal
    tumor_filtered = get_filtered_per_alias(tumor, tumor_alias)
    normal_filtered = get_filtered_per_alias(normal, "normal")
    all_filtered = (
        pd.concat([tumor_filtered, normal_filtered])
        .groupby("id", group_keys=False)
        .apply(diff_tumor_normal_peptides, column="pep_seq", tumor_alias=tumor_alias)
    )
    info_tidy = tidy_info(info, tumor_alias)

    # netMHCpan 4.1 truncates the fasta entry IDs, so we have to cut down the IDs
    # that microphaser originally provided to make the following .join() work
    len_tumor_id = len(tumor_filtered["id"][0])
    len_normal_id = len(normal_filtered["id"][0])
    assert len_tumor_id == len_normal_id, f"'id's' are of different length, tumor: {len_tumor_id}, normal: {len_normal_id}, please check your input data.\n"
    info_tidy["id"] = info_tidy["id"].str[:len_tumor_id]
    # Double-check for duplicates resulting from the id truncation
    check_duplicates(info_tidy, ["id", "alias"], specific_error="Here, the problem is most likely the truncation of 'id's by netMHCpan.\n")

    all_annotated = all_filtered.merge(info_tidy, how="left", on=["id", "alias"])

    # Double-check for weird duplicates, as previously done in Jan's code.
    check_duplicates(all_annotated, ["transcript", "offset", "pep_seq", "aa_changes"])

    column_order = [
        "id",
        "pep_seq",
        "pos_in_id_seq",
        "alias",
        "num_binders",
        "freq",
        "depth",
        "num_var_sites",
        "num_var_in_pep",
        "top_el_rank_allele",
        "top_el_rank_bind_core",
        "top_el_rank_el_rank",
        "top_el_rank_el_score",
        "ave_el_score",
        "top_ba_rank_allele",
        "top_ba_rank_bind_core",
        "top_ba_rank_ba_rank",
        "top_ba_rank_ba_score",
        "aa_changes",
        "genomic_pos",
        "nt_seq",
        "gene_name",
        "gene_id",
        "transcript",
        "chrom",
        "offset",
        "frame",
        "strand",
        "freq_credible_interval",
    ]

    def get_id_rank(group: pd.DataFrame, tumor_alias: str):
        return group.loc[group["alias"] == tumor_alias, 'top_el_rank_el_rank'].squeeze()
    
    sort_rank = all_annotated.groupby('id').apply(get_id_rank, tumor_alias).rename('sort_rank')
    all_sorted = all_annotated.merge(sort_rank, on=['id'], how='left').sort_values(['sort_rank', 'id', 'alias'], ascending=[True, True, False]).drop(columns='sort_rank')

    return all_sorted.reindex(columns=column_order)


info = pd.read_csv(snakemake.input.info, sep="\t", dtype=str)
tumor = pd.read_csv(snakemake.input.neo, sep="\t", dtype={"pos_in_id_seq": str})
normal = pd.read_csv(snakemake.input.normal, sep="\t", dtype={"pos_in_id_seq": str})
data = merge_data_frames(info, tumor, normal, snakemake.wildcards.tumor_alias)
data.to_csv(snakemake.output[0], index=False, sep="\t")
