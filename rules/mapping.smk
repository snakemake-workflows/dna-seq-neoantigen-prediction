def get_reads(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("../raw/{sample}_{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "raw/{sample}.fastq.gz".format(**wildcards)

rule map:
    input:
        reads=get_reads
    output:
        "bwa/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
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
        "../bwa/{sample}.bam"
    output:
        bam="bwa/{sample}.rmdup.bam",
        metrics="dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.31.1/bio/picard/markduplicates"

rule index:
    input:
        bam="bwa/{sample}.rmdup.bam"
    output:
        bam="bwa/{sample}.rmdup.bam.bai"
    params: ""
    wrapper:
        "0.31.1/bio/samtools/index"
