rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    log:
        "logs/bcf-index/{prefix}.log"
    wrapper:
        "0.60.0/bio/bcftools/index"        


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.59.2/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz"
    output:
        "resources/{prefix}.{format}.gz.tbi"
    log:
        "logs/tabix/{prefix}.{format}.log"
    params:
        get_tabix_params
    cache: True
    wrapper:
        "0.59.2/bio/tabix"

rule gzip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fastq.gz"
    log:
        "logs/gz-fastq/{prefix}.log"
    shell:
        "gzip < {input} > {output}"