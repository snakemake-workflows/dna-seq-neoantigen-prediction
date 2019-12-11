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
    return [samples.loc[(wildcards.group), "matched_normal"], samples.loc[wildcards.group, "sample"]]

def get_paired_bams(wildcards):
    return expand("bwa/{sample}.rmdup.bam", sample=get_paired_samples(wildcards))

def get_paired_bais(wildcards):
    return expand("bwa/{sample}.rmdup.bam.bai", sample=get_paired_samples(wildcards))

def get_normal(wildcards):
    return samples.loc[(wildcards.sample), "matched_normal"]

def get_reads(wildcards):
    return get_seperate(wildcards.sample, wildcards.group)

def get_seperate(sample, group):
    return units.loc[(sample, "DNA"), "fq" + str(group)]

def get_proteome(wildcards):
    return expand("microphaser/fasta/germline/{normal}/reference_proteome.bin", normal=get_normal(wildcards))#samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_optitype(wildcards):
    return expand("optitype/{germline}/hla_alleles_{germline}.tsv", germline=get_normal(wildcards))

def get_germline_hla(wildcards):
    return expand("HLA-LA/hlaII_{germline}.tsv", germline=get_normal(wildcards))

def get_normal_bam(wildcards):
    return expand("bwa/{normal}.rmdup.bam", normal=get_normal(wildcards))

def get_normal_bai(wildcards):
    return expand("bwa/{normal}.rmdup.bam.bai", normal=get_normal(wildcards))

def get_germline_variants(wildcards):
    return expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz", germline=get_normal(wildcards))

def get_germline_variants_index(wildcards):
    return expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz.csi", germline=get_normal(wildcards))

def get_group_observations(wildcards):
    return expand("observations/{group}/{sample}.{caller}.{contig}.bcf", 
                  contig=wildcards.contig, 
                  caller=wildcards.caller, 
                  group=wildcards.group,
                  sample=get_paired_samples(wildcards))

def get_group_aliases(wildcards):
    return samples.loc[(samples["sample"] == wildcards.group) | (samples["sample"] == samples.loc[wildcards.group, "matched_normal"])]["type"]

wildcard_constraints:
    group="|".join(samples[samples.type == "tumor"]["sample"]),
    sample="|".join(samples["sample"]),
    contig="|".join(contigs),
    caller="|".join(["freebayes", "delly"])

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c["activate"])

caller=list(filter(None, ["freebayes" if is_activated("calling/freebayes") else None, "delly" if is_activated("calling/delly") else None]))
