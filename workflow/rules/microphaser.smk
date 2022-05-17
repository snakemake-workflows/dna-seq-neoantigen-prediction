rule norm_bcf:
    input:
        "results/final-calls/{group}.{set}.bcf",
        genome="resources/genome.fasta",
    output:
        "results/final-calls/{group}.{set}.norm.bcf",
    log:
        "logs/bcftools/norm/{group}.{set}.log",
    params:
        lambda w, input: "-f {} -O b -m-".format(input.genome),  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.65.0/bio/bcftools/norm"


rule merge_tumor_normal:
    input:
        calls=expand(
            "results/final-calls/{{group}}.{sets}.norm.bcf",
            sets=[
                config["params"]["microphaser"]["variant_sets"]["normal"],
                config["params"]["microphaser"]["variant_sets"]["tumor"],
            ],
        ),
        index=expand(
            "results/final-calls/{{group}}.{sets}.norm.csi",
            sets=[
                config["params"]["microphaser"]["variant_sets"]["normal"],
                config["params"]["microphaser"]["variant_sets"]["tumor"],
            ],
        ),
    output:
        "results/final-calls/{group}.merged_tumor_normal.norm.bcf",
    log:
        "bcftools/concat-tumor-normal/{group}.merged_tumor_normal.log",
    params:
        extra="-O b -a",
    wrapper:
        "0.64.0/bio/bcftools/concat"


rule microphaser_tumor:
    input:
        bcf="results/final-calls/{group}.merged_tumor_normal.norm.annotated.bcf",
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        mt_fasta="results/microphaser/fasta/{group}/{tumor_alias}.merged_tumor_normal.{contig}.neo.fa",
        wt_fasta="results/microphaser/fasta/{group}/{tumor_alias}.merged_tumor_normal.{contig}.normal.fa",
        tsv="results/microphaser/info/{group}/{tumor_alias}.merged_tumor_normal.{contig}.tsv",
    log:
        "logs/microphaser_tumor/{group}/{tumor_alias}.merged_tumor_normal.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser somatic {input.bam} --variants {input.bcf} --ref {input.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.window_length} "
        "< {input.track} > {output.mt_fasta} 2> {log}"


rule microphaser_normal:
    input:
        bcf="results/final-calls/{group}.{normal_set}.variants.reheader.norm.bcf",
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        wt_fasta=("results/microphaser/fasta/{group}/{normal_alias}.{normal_set}.{contig}.fa"),
        wt_tsv=("results/microphaser/info/{group}/{normal_alias}.{normal_set}.{contig}.tsv"),
    log:
        "logs/microphaser_germline/{group}/{normal_alias}.{normal_set}-{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser normal {input.bam} --variants {input.bcf} --ref {input.ref} -t {output.wt_tsv} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta} 2> {log}"


rule concat_normal_proteome:
    input:
        expand(
            "results/microphaser/fasta/{{group}}/normal.{{normal_set}}.{contig}.fa",
            contig=contigs,
        ),
    output:
        "results/microphaser/fasta/{group}.{normal_set}.normal_proteome.fa",
    log:
        "logs/microphaser/concat_normal_proteome/{group}.{normal_set}.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule build_normal_proteome_db:
    input:
        "results/microphaser/fasta/{group}.{normal_set}.normal_proteome.fa",
    output:
        bin="results/microphaser/bin/{group}.{normal_set}.{mhc}.normal_proteome.bin",
        fasta="results/microphaser/fasta/{group}.{normal_set}.{mhc}.normal_proteome.peptides.fasta",
    log:
        "logs/microphaser/build_normal_proteome_db/{group}.{normal_set}-{mhc}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        length=lambda wildcards: config["params"]["microphaser"]["peptide_len"][
            wildcards.mhc
        ],
    shell:
        "microphaser build_reference -r {input} -o {output.bin} -l {params.length} --peptides {output.fasta} > {log} 2>&1"


rule microphaser_filter:
    input:
        tsv="results/microphaser/info/{group}/{tumor_alias}.merged_tumor_normal.{contig}.tsv",
        proteome=expand(
            "results/microphaser/bin/{{group}}.{normal_set}.{{mhc}}.normal_proteome.bin",
            normal_set=config["params"]["microphaser"]["variant_sets"]["normal"],
        ),
    output:
        mt_fasta=(
            "results/microphaser/fasta/filtered/{group}/{tumor_alias}.merged_tumor_normal.{mhc}.{contig}.neo.fa"
        ),
        wt_fasta=(
            "results/microphaser/fasta/filtered/{group}/{tumor_alias}.merged_tumor_normal.{mhc}.{contig}.normal.fa"
        ),
        tsv="results/microphaser/info/filtered/{group}/{tumor_alias}.merged_tumor_normal.{mhc}.{contig}.tsv",
        removed="results/microphaser/info/removed/{group}/{tumor_alias}.merged_tumor_normal.{mhc}.{contig}.removed.tsv",
    log:
        "logs/microphaser_filter/{group}/{tumor_alias}.merged_tumor_normal.{mhc}.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        length=lambda wildcards: config["params"]["microphaser"]["peptide_len"][
            wildcards.mhc
        ],
    shell:
        "microphaser filter -r {input.proteome} -t {input.tsv} -o {output.tsv} -n {output.wt_fasta} -s {output.removed} -l {params.length} > {output.mt_fasta} 2>{log}"


rule concat_tsvs:
    input:
        expand(
            "results/microphaser/info/filtered/{{group}}/{{tumor_alias}}.merged_tumor_normal.{{mhc}}.{contig}.tsv",
            contig=contigs,
        ),
    output:
        "results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.tsv",
    log:
        "logs/concat_tsvs/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output} 2>{log}"
