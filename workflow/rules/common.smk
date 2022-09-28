import glob

import pandas as pd
from os import path
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

##### config file #####

validate(config, schema="../schemas/config.schema.yaml")

##### samples sheet #####

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={
            "sample_name": str,
            "group": str,
            "alias": str,
        },
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

##### units sheet #####

units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={
            "sample_name": str,
            "sequencing_type": str,
            "unit_name": str,
            "adapters": str,
        },
        comment="#",
    )
    .set_index(["sample_name", "sequencing_type", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

groups = samples["group"].unique()

if "groups" in config:
    group_annotation = (
        pd.read_csv(config["groups"], sep="\t", dtype={"group": str})
        .set_index("group")
        .sort_index()
    )
    group_annotation = group_annotation.loc[groups]
else:
    group_annotation = pd.DataFrame({"group": groups}).set_index("group")


contigs = [c for c in range(1, 23)]
contigs.extend(["X", "Y"])

##### constraining wildcards #####


wildcard_constraints:
    cancer_sample="|".join(samples[samples.alias != "normal"]["sample_name"]),
    normal_sample="|".join(samples[samples.alias == "normal"]["sample_name"]),
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),
    alias="|".join(pd.unique(samples["alias"])),
    tumor_alias="|".join(
        pd.unique(samples.loc[samples["alias"].str.match("tumor"), "alias"])
    ),
    normal_alias="normal",
    tumor_set=config["params"]["microphaser"]["events"]["tumor"],
    normal_set=config["params"]["microphaser"]["events"]["normal"],
    set="|".join(
        [
            config["params"]["microphaser"]["events"]["tumor"],
            config["params"]["microphaser"]["events"]["normal"],
        ]
    ),
    group="|".join(pd.unique(samples["group"])),
    caller="|".join(["freebayes", "delly"]),
    peptide_type="|".join(["normal", "neo"]),
    event="|".join(["somatic", "germline", "complete"]),
    read="|".join(["single", "R1", "R2"]),
    seqtype="|".join(["DNA", "RNA"]),


### Output generation ###


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_final_output():
    final_output = []
    for group in pd.unique(samples["group"]):
        smps = samples.loc[samples["group"] == group, "sample_name"]
        tumor_aliases = samples.loc[
            (samples["group"] == group) & (samples["alias"].str.match("tumor")),
            "alias",
        ]
        if config["neoantigen_prediction"]["activate"]:
            final_output.extend(
                expand(
                    "results/datavzrd/neoprint/{group}.{tumor_alias}.{mhc}",
                    group=group,
                    tumor_alias=tumor_aliases,
                    mhc=["I", "II"],
                )
            )
        #    sequencing_types = pd.unique(
        #       units.loc[units["sample_name"].isin(smps), "sequencing_type"]
        #    )
        #    final_output.extend(
        #       expand(
        #           "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.{seqtype}.tsv",
        #           group=group,
        #           tumor_alias=tumor_aliases,
        #           mhc=list(
        #               filter(
        #                   None,
        #                   [
        #                       "net_mhc_pan"
        #                       if is_activated("params/net_mhc_pan")
        #                       else None,
        #                       "net_mhc_two_pan"
        #                       if is_activated("params/net_mhc_two_pan")
        #                       else None,
        #                   ],
        #               )
        #           ),
        #           seqtype=sequencing_types,
        #       )
        #    )
        else:
            final_output = expand(
                [
                    "results/hla_la/{group}.{tumor_alias}.hlaI.tsv",
                    "results/hla_la/{group}.{tumor_alias}.hlaII.tsv",
                ],
                group=group,
                tumor_alias=tumor_aliases,
            )

    return final_output


caller = list(
    filter(None, ["freebayes" if is_activated("calling/freebayes") else None])
)

### helper functions ###


def is_paired_end(sample, seqtype):
    sample_units = units.loc[sample].loc[seqtype]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {} and sequencing type {}, must be all paired end or all single end".format(
        sample, seqtype
    )
    return all_paired


def get_sample_from_group_and_alias(group, alias):
    sample = samples.loc[
        (samples["group"] == group) & (samples["alias"] == alias), "sample_name"
    ].squeeze()
    return sample


def get_bam_from_group_and_alias(ext=".bam"):
    def inner(wildcards):
        alias = wildcards.get(
            "alias",
            wildcards.get("tumor_alias", wildcards.get("normal_alias", "unknown")),
        )
        if alias == "unknown":
            raise CustomException(
                "get_bam_from_group_and_alias() requires one of the following wildcards: 'alias', 'tumor_alias', 'normal_alias'."
            )
        sample = get_sample_from_group_and_alias(wildcards.group, alias)
        return f"results/recal/{sample}.sorted{ext}"

    return inner


def get_alleles_MHCI(wildcards):
    alias = "normal" if wildcards.peptide_type == "normal" else wildcards.tumor_alias
    return expand(
        "results/hla_la/{group}.{alias}.hlaI.tsv",
        group=wildcards.group,
        alias=alias,
    )


def get_alleles_MHCII(wildcards):
    alias = "normal" if wildcards.peptide_type == "normal" else wildcards.tumor_alias
    return expand(
        # TODO: check that hlaII is correct here, and not hlaI which it previously was
        "results/hla_la/{group}.{alias}.hlaII.tsv",
        group=wildcards.group,
        alias=alias,
    )


##### Other stuff ####

neofox_important_cols = {
    "general": [
            "gene",
            "mutation_mutatedXmer",
            "mutation_wildTypeXmer",
            "purity_adjusted_DNA_VAF",
            "imputedGeneExpression",
        ],
    "I": [
            "PRIME_best_rank",
            "PRIME_best_score",
            "PRIME_best_peptide",
            "PRIME_best_allele",
            "Recognition_Potential_MHCI_9mer",
            "Improved_Binder_MHCI",
            "Selfsimilarity_MHCI_conserved_binder",
            "Best_rank_MHCI_9mer_score",
            "Best_rank_MHCI_9mer_score_WT",
            "Best_rank_MHCI_9mer_epitope",
            "Best_rank_MHCI_9mer_epitope_WT",
            "Best_rank_MHCI_9mer_allele",
            "Best_rank_MHCI_9mer_allele_WT",
        ],
    "II": [
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
}