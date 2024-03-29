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
