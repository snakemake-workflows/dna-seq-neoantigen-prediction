def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
#    if not pd.isnull(units.loc[(wildcards.sample), "rna_fq1"]):
#        return units.loc[
#            (wildcards.sample), ["rna_fq1", "rna_fq2"]].dropna()
    return get(wildcards.sample, "RNA")

rule kallisto_index:
    input:
        "resources/genome.cdna.fasta"
    output:
        "resources/kallisto/transcripts.idx"
    params:
        extra=""
    log:
        "results/logs/kallisto/index.log"
    cache: True
    wrapper:
        "0.60.1/bio/kallisto/index"


def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if len(input.fq) == 1:
        extra += " --single"
        extra += (" --fragment-length {unit.fragment_len_mean} "
                  "--sd {unit.fragment_len_sd}").format(
                    unit=units.loc[
                        (wildcards.sample, wildcards.unit)])
    else:
        extra += " --fusion"
    return extra


rule kallisto_quant:
    input:
        fastq=get_fastqs,
        index="resources/kallisto/transcripts.idx"
    output:
        directory("results/kallisto/{sample}")
    params:
        extra=kallisto_params
    log:
        "results/logs/kallisto/quant/{sample}.log"
    wrapper:
        "0.60.1/bio/kallisto/quant"
