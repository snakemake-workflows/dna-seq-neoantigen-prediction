import glob

import pandas as pd
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
    if config["epitope_prediction"]["activate"]:
        for group in pd.unique(samples["group"]):
            smps = samples.loc[samples["group"] == group, "sample_name"]
            sequencing_types = pd.unique(
                units.loc[units["sample_name"].isin(smps), "sequencing_type"]
            )
            tumor_aliases = samples.loc[
                (samples["group"] == group) & (samples["alias"].str.match("tumor")),
                "alias",
            ]
            final_output.extend(
                expand(
                    "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.{seqtype}.tsv",
                    group=group,
                    tumor_alias=tumor_aliases,
                    mhc=list(
                        filter(
                            None,
                            [
                                "netMHCpan"
                                if is_activated("affinity/netMHCpan")
                                else None,
                                "netMHCIIpan"
                                if is_activated("affinity/netMHCIIpan")
                                else None,
                            ],
                        )
                    ),
                    seqtype=sequencing_types,
                )
            )
    else:
        if config["HLAtyping"]["HLA_LA"]["activate"]:
            final_output = expand(
                [
                    "results/optitype/{group}/{group}.{alias}.hla_alleles.tsv",
                    "results/HLA-LA/{group}.{alias}.hlaI.tsv",
                    "results/HLA-LA/{group}.{alias}.hlaII.tsv",
                ],
                sample=samples["sample_name"],
            )
        else:
            final_output = expand(
                "results/optitype/{group}/{group}.{alias}.hla_alleles.tsv",
                sample=samples["sample_name"],
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


def get_optitype_reads_input(wildcards):
    sample = get_sample_from_group_and_alias(wildcards.group, wildcards.alias)
    if is_activated("HLAtyping/optitype_prefiltering"):
        if is_paired_end(sample, "DNA"):
            return expand(
                "results/razers3/fastq/{sample}_{read}.fished.fastq",
                sample=sample,
                read=["R1", "R2"],
            )
        return f"results/razers3/fastq/{sample}_single.fastq"
    else:
        wildcards["sample"] = sample
        return get_map_reads_input(wildcards)


def get_bam_from_group_and_alias(ext=".bam"):
    def inner(wildcards):
        alias = wildcards.get("alias",
            wildcards.get("tumor_alias",
                wildcards.get("normal_alias", "unknown")
            )
        )
        if alias == "unknown":
            raise CustomException(
                "get_bam_from_group_and_alias() requires on of the following wildcards: 'alias', 'tumor_alias', 'normal_alias'."
            )
        sample = get_sample_from_group_and_alias(wildcards.group, alias)
        return f"results/recal/{sample}.sorted{ext}"

    return inner


def get_alleles_MHCI(wildcards):
    alias = "normal" if wildcards.peptide_type == "normal" else wildcards.tumor_alias
    return expand(
        "results/optitype/{group}/{group}.{alias}.hla_alleles.tsv",
        group=wildcards.group,
        alias=alias,
    )


def get_alleles_MHCII(wildcards):
    alias = "normal" if wildcards.peptide_type == "normal" else wildcards.tumor_alias
    return expand(
        #TODO: check that hlaII is correct here, and not hlaI which it previously was
        "results/HLA-LA/{group}.{alias}.hlaII.tsv", group=wildcards.group, alias=alias
    )
