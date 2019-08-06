def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
#    if not pd.isnull(units.loc[(wildcards.sample), "rna_fq1"]):
#        return units.loc[
#            (wildcards.sample), ["rna_fq1", "rna_fq2"]].dropna()
    return get(wildcards.sample, "RNA")

rule kallisto_index:
    input:
        config["reference"]["transcriptome"]
    output:
        "transcriptome/kallisto/transcripts.idx"
    log:
        "logs/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"


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
        fq=get_fastqs,
        idx="transcriptome/kallisto/transcripts.idx"
    output:
        abundance="transcriptome/kallisto/{sample}/abundance.tsv"
    log:
        "logs/kallisto/quant/{sample}.log"
    params:
        extra=kallisto_params
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o transcriptome/kallisto/{wildcards.sample} "
        "{params.extra} {input.fq} 2> {log}"
