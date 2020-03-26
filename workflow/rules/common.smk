import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### sample sheets #####

samples = pd.read_csv(config["samples"], sep='\t').set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep='\t').set_index(["sample", "sequencing_type"], drop=False)
print(units)
#validate(units, schema="schemas/units.schema.yaml")

contigs = pd.read_csv(
    config["reference"]["genome"] + ".fai",
    header=None, usecols=[0], squeeze=True, dtype=str, sep="\t"
)

# Use this to ignore decoy and unplaced contigs.
contigs = contigs[contigs.str.contains("_|M") == False]

def get_DNA_reads(wildcards):
    if config["trimming"]["skip"]:
        # no need for trimming, return raw reads
        return units.loc[(wildcards.sample, "DNA"), ["fq1", "fq2"]]

def get_paired_samples(wildcards):
    return [samples.loc[(wildcards.pair), "matched_normal"], samples.loc[wildcards.pair, "sample"]]

def get_paired_bams(wildcards):
    return expand("results/bwa/{sample}.rmdup.bam", sample=get_paired_samples(wildcards))

def get_paired_bais(wildcards):
    return expand("results/bwa/{sample}.rmdup.bam.bai", sample=get_paired_samples(wildcards))

def get_normal(wildcards):
    return samples.loc[(wildcards.sample), "matched_normal"]

def get_reads(wildcards):
    return get_seperate(wildcards.sample, wildcards.group)

def get_seperate(sample, group):
    return units.loc[(sample, "DNA"), "fq" + str(group)]

def get_proteome(wildcards):
    return expand("results/microphaser/fasta/germline/{normal}/{mhc}/reference_proteome.bin", normal=get_normal(wildcards), mhc=wildcards.mhc)

def get_germline_optitype(wildcards):
    return expand("results/optitype/{germline}/hla_alleles_{germline}.tsv", germline=get_normal(wildcards))

def get_germline_hla(wildcards):
    return expand("results/HLA-LA/hlaII_{germline}.tsv", germline=get_normal(wildcards))

def get_normal_bam(wildcards):
    return expand("results/bwa/{normal}.rmdup.bam", normal=get_normal(wildcards))

def get_normal_bai(wildcards):
    return expand("results/bwa/{normal}.rmdup.bam.bai", normal=get_normal(wildcards))

def get_germline_variants(wildcards):
    return expand("results/strelka/germline/{germline}/results/variants/variants.reheader.bcf", germline=get_normal(wildcards))

def get_germline_variants_index(wildcards):
    return expand("results/strelka/germline/{germline}/results/variants/variants.reheader.bcf.csi", germline=get_normal(wildcards))

def get_pair_observations(wildcards):
    return expand("results/observations/{pair}/{sample}.{caller}.bcf", 
                  caller=wildcards.caller, 
                  pair=wildcards.pair,
                  sample=get_paired_samples(wildcards))

def get_pair_aliases(wildcards):
    return samples.loc[(samples["sample"] == wildcards.pair) | (samples["sample"] == samples.loc[wildcards.pair, "matched_normal"])]["type"]

wildcard_constraints:
    pair="|".join(samples[samples.type == "tumor"]["sample"]),
    sample="|".join(samples["sample"]),
    caller="|".join(["freebayes", "delly"])

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c["activate"])

caller=list(filter(None, ["freebayes" if is_activated("calling/freebayes") else None, "delly" if is_activated("calling/delly") else None]))
