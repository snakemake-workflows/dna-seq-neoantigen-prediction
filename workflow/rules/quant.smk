def kallisto_params(wildcards, input):
    extra = config["params"]["kallisto"]
    if len(input.fastq) == 1:
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
        fastq=get_quant_reads_input,
        index="resources/kallisto/transcripts.idx"
    output:
        directory("results/kallisto/{sample}")
    params:
        extra=kallisto_params
    log:
        "results/logs/kallisto/quant/{sample}.log"
    wrapper:
        "0.60.1/bio/kallisto/quant"
