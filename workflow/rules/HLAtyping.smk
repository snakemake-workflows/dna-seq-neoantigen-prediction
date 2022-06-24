rule HLA_LA:
    input:
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        index="resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
    output:
        "results/HLA-LA/output/{group}_{alias}/hla/R1_bestguess_G.txt",
    threads: 7
    log:
        "logs/HLA-LA/{group}_{alias}.log",
    params:
        graph=lambda w, input: os.path.basename(os.path.dirname(input.index)),
        graphdir=lambda w, input: os.path.dirname(os.path.dirname(input.index)),
    conda:
        "../envs/hla_la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.group}_{wildcards.alias} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir results/HLA-LA/output --maxThreads {threads} > {log} 2>&1"


rule parse_HLA_LA:
    input:
        "results/HLA-LA/output/{group}_{alias}/hla/R1_bestguess_G.txt",
    output:
        report(
            "results/HLA-LA/{group}.{alias}.hlaI.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
        report(
            "results/HLA-LA/{group}.{alias}.hlaII.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
    log:
        "logs/parse-HLA-LA/{group}.{alias}.log",
    script:
        "../scripts/parse_HLA_types.py"


rule get_hla_aligning_reads:
    input:
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        regions="resources/hla_alleles/hla_allele_regions.expanded_1000.bed",
    output:
        bam="results/fished/{group}.{alias}.bam",
        idx="results/fished/{group}.{alias}.bai",
    log:
        "logs/get_hla_reads/{group}.{alias}.log",
    params:
        extra=lambda wc, input: f"--regions-file {input.regions}"
    wrapper:
        "v1.7.0/bio/samtools/view"


ruleorder: get_hla_aligning_reads > bam_index


rule hla_reads_single_ends:
    input:
        "results/fished/{group}.{alias}.bam",
        "results/fished/{group}.{alias}.bai",
    output:
        bam="results/fished/{group}.{alias}.{read}.bam",
        idx="results/fished/{group}.{alias}.{read}.bai",
    log:
        "logs/split_hla_reads/{group}.{alias}.{read}.log",
    params:
        extra=lambda wc: "-f 0x80" if wc.read == "R2" else "-f 0x40"
    wrapper:
        "v1.7.0/bio/samtools/view"


ruleorder: hla_reads_single_ends > bam_index
