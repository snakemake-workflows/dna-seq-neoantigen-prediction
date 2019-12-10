def get_DNA_reads(wildcards):
    if config["trimming"]["skip"]:
        # no need for trimming, return raw reads
        return units.loc[(wildcards.sample, "DNA"), ["fq1", "fq2"]]

def get_reads(wildcards):
    return get_seperate(wildcards.sample, wildcards.group)

def get_seperate(sample, group):
    return units.loc[(sample, "DNA"), "fq" + str(group)]

def get_proteome(wildcards):
    return(expand("microphaser/fasta/germline/{normal}/reference_proteome.bin", normal=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_optitype(wildcards):
    return(expand("optitype/{germline}/hla_alleles_{germline}.tsv", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_hla(wildcards):
    return(expand("HLA-LA/hlaII_{germline}.tsv", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_normal_bam(wildcards):
    return(expand("bwa/{normal}.rmdup.bam", normal=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_normal_bai(wildcards):
    return(expand("bwa/{normal}.rmdup.bam.bai", normal=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_variants(wildcards):
    return(expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_variants_index(wildcards):
    return(expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz.csi", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))
