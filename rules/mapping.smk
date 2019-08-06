def get_reads(wildcards):
    if config["trimming"]["skip"]:
        # no need for trimming, return raw reads
        return units.loc[(wildcards.sample, wildcards.typ), ["fq1", "fq2"]].tolist()
    else:
        # trimming needed, use trimmed reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("trimmed/{sample}-{typ}.fastq.gz",
                    group=[1,2], **wildcards)
        # single end sample
        return "trimmed/{sample}-{typ}.fastq.gz".format(**wildcards)

rule map:
    input:
        reads=get_reads
    output:
        "bwa/{sample}-{typ}.bam"
    log:
        "logs/bwa_mem/{sample}-{typ}.log"
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
        "bwa/{sample}-{typ}.bam"
    output:
        bam="bwa/{sample}-{typ}.rmdup.bam",
        metrics="dedup/{sample}-{typ}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{typ}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.31.1/bio/picard/markduplicates"

rule index:
    input:
        bam="bwa/{sample}-{typ}.rmdup.bam"
    output:
        bam="bwa/{sample}-{typ}.rmdup.bam.bai"
    params: ""
    wrapper:
        "0.31.1/bio/samtools/index"
