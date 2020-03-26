rule map:
    input:
        reads=get_DNA_reads
    output:
        "results/bwa/{sample}.bam"
    log:
        "results/logs/bwa_mem/{sample}.log"
    params:
        index=config["reference"]["genome"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",
        sort_order="coordinate",
        sort_extra=""
    threads: 8
    wrapper:
        "0.31.1/bio/bwa/mem"

rule mark_duplicates:
    input:
        "results/bwa/{sample}.bam"
    output:
        bam="results/bwa/{sample}.rmdup.bam",
        metrics="results/dedup/{sample}.metrics.txt"
    log:
        "results/logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.31.1/bio/picard/markduplicates"

rule index:
    input:
        bam="results/bwa/{sample}.rmdup.bam"
    output:
        bam="results/bwa/{sample}.rmdup.bam.bai"
    params: ""
    wrapper:
        "0.31.1/bio/samtools/index"
