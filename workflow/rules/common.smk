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
                    "results/neoantigens/{group}.{tumor_alias}.{tumor_event}.{mhc}.{seqtype}.tsv",
                    group=group,
                    tumor_alias=tumor_aliases,
                    tumor_event=config["params"]["microphaser"]["events"]["tumor"],
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
                    "results/optitype/{sample}/hla_alleles_{sample}.tsv",
                    "results/HLA-LA/hlaI_{sample}.tsv",
                    "results/HLA-LA/hlaII_{sample}.tsv",
                ],
                sample=samples["sample_name"],
            )
        else:
            final_output = expand(
                "results/optitype/{sample}/hla_alleles_{sample}.tsv",
                sample=samples["sample_name"],
            )
    return final_output


def get_fusion_output():
    if config["fusion"]["arriba"]["activate"]:
        fusion_output = expand(
            "results/fusion/arriba/{sample}.fusions.tsv",
            sample=units[units["sequencing_type"] == "RNA"]["sample_name"],
        )
    else:
        fusion_output = []
    return fusion_output


def get_tmb_targets():
    if is_activated("tmb"):
        return expand(
            "results/plots/tmb/{group}.{mode}.svg",
            group=samples[(samples.alias == "tumor")]["sample_name"],
            mode=config["tmb"].get("mode", "curve"),
        )
    else:
        return []


caller = list(
    filter(None, ["freebayes" if is_activated("calling/freebayes") else None])
)

### helper functions ###

## alignment ##


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.seqtype].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{T}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name,
            U=unit.unit_name,
            T=unit.sequencing_type,
            E=ending,
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{T}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name,
                U=unit.unit_name,
                T=unit.sequencing_type,
                E=ending,
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    pattern = (
        units.loc[wildcards.sample]
        .loc[wildcards.seqtype]
        .loc[wildcards.unit, wildcards.fq]
    )
    if "*" in pattern:
        files = sorted(
            glob.glob(
                units.loc[wildcards.sample]
                .loc[wildcards.seqtype]
                .loc[wildcards.unit, wildcards.fq]
            )
        )
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}, sequencing type {}). "
                "Please check the your sample sheet.".format(
                    wildcards.unit, wildcards.sample, wildcards.seqtype
                )
            )
    else:
        files = [pattern]
    return files


def get_cutadapt_adapters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        adapters = unit["adapters"]
        if isinstance(adapters, str):
            return adapters
        return ""
    except KeyError:
        return ""


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


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}/{seqtype}/{unit}_{read}.fastq.gz",
            unit=units.loc[
                (units["sequencing_type"] == wc.seqtype)
                & (units["sample_name"] == wc.sample),
                "unit_name",
            ],
            sample=wc.sample,
            read=wc.read,
            seqtype=wc.seqtype,
        )
    unit = units.loc[wc.sample].loc[wc.seqtype]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand(
            "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return units.loc[wc.sample].loc[wc.seqtype, fq].tolist()


def get_map_reads_input(sample):
    if is_paired_end(sample, "DNA"):
        return [
            "results/merged/DNA/{sample}_R1.fastq.gz",
            "results/merged/DNA/{sample}_R2.fastq.gz",
        ]
    return "results/merged/DNA/{sample}_single.fastq.gz"


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, platform=samples.loc[wildcards.sample, "platform"]
    )


def get_recalibrate_quality_input(bai=False):
    ext = ".bai" if bai else ""

    def inner(wildcards):
        if is_activated("remove_duplicates"):
            return f"results/dedup/{wildcards.sample}.sorted.bam{ext}"
        else:
            return f"results/mapped/{wildcards.sample}.sorted.bam{ext}"

    return inner


def get_optitype_reads_input(wildcards):
    sample = get_sample_from_group_and_alias(wildcards.group, wildcards.alias)
    if is_activated("HLAtyping/optitype_prefiltering"):
        if is_paired_end(sample, "DNA"):
            return expand(
                "results/razers3/fastq/{sample}_{read}.fished.fastq",
                sample=sample,
                read=["R1", "R2"],
            )
        return "results/razers3/fastq/{sample}_single.fastq"
    else:
        return get_map_reads_input(sample)


def get_oncoprint_batch(wildcards):
    if wildcards.batch == "all":
        groups = samples[samples["alias"] == "tumor"]["sample_name"].unique()
    else:
        groups = samples.loc[
            samples[config["oncoprint"]["stratify"]["by-column"]] == wildcards.batch,
            "group",
        ].unique()
    return expand(
        "results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=groups
    )


def get_tabix_params(wildcards):
    if wildcards.format == "vcf":
        return "-p vcf"
    if wildcards.format == "txt":
        return "-s 1 -b 2 -e 2"
    raise ValueError("Invalid format for tabix: {}".format(wildcards.format))


def get_sample_from_group_and_alias(group, alias):
    sample = samples.loc[
        (samples["group"] == group) & (samples["alias"] == alias), "sample_name"
    ].squeeze()
    return sample


def get_normal_bam(ext=".bam"):
    def inner(wildcards):
        normal_sample = get_sample_from_group_and_alias(wildcards.group, "normal")
        return f"results/recal/{normal_sample}.sorted{ext}"

    return inner


def get_tumor_bam_from_group_and_alias(ext=".bam"):
    def inner(wildcards):
        tumor_sample = get_sample_from_group_and_alias(
            wildcards.group, wildcards.tumor_alias
        )
        return f"results/recal/{tumor_sample}.sorted{ext}"

    return inner


def get_bam_from_group_and_alias(ext=".bam"):
    def inner(wildcards):
        sample = get_sample_from_group_and_alias(wildcards.group, wildcards.alias)
        return f"results/recal/{sample}.sorted{ext}"

    return inner


def get_quant_reads_input(wildcards):
    sample = get_sample_from_group_and_alias(wildcards.group, wildcards.tumor_alias)
    if is_paired_end(sample, "RNA"):
        return [
            "results/merged/RNA/{sample}_R1.fastq.gz",
            "results/merged/RNA/{sample}_R2.fastq.gz",
        ]
    return "results/merged/RNA/{sample}_single.fastq.gz"


def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if len(input.fastq) == 1:
        extra += " --single"
        extra += (
            " --fragment-length {unit.fragment_len_mean} " "--sd {unit.fragment_len_sd}"
        ).format(unit=units.loc[(wildcards.sample, wildcards.unit)])
    else:
        extra += " --fusion"
    return extra


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
        "results/HLA-LA/{group}.{alias}.hlaI.tsv", group=wildcards.group, alias=alias
    )
