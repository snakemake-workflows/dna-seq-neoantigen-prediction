rule samtools_index:
    input:
        bam="{bam}.bam"
    output:
        bam="{bam}.bam.bai"
    params: ""
    wrapper:
        "0.31.1/bio/samtools/index"

rule map_reads:
    input:
        reads=get_DNA_reads,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="resources/genome.fasta",
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.39.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.39.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam="results/dedup/{sample}.sorted.bam",
        bai="results/dedup/{sample}.sorted.bam.bai",
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=("results/recal/{sample}.grp")
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"]
        java_opts=""
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    threads: 8
    wrapper:
        "0.62.0/bio/gatk/baserecalibrator"

rule apply_bqsr:
    input:
        bam="results/dedup/{sample}.sorted.bam",
        bai="results/dedup/{sample}.sorted.bam.bai",
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}.grp"
    output:
        bam="results/recal/{sample}.sorted.bam"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    wrapper:
        "0.62.0/bio/gatk/applybqsr"