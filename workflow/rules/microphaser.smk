rule microphaser_tumor:
    input:
        vcf="results/strelka/merged/{group}.{tumor_event}.{normal_event}.norm.annotated.bcf",
        bam=get_tumor_bam(),
        bai=get_tumor_bam(ext=".bam.bai"),
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        mt_fasta="results/microphaser/fasta/{group}/tumor.{tumor_event}.{normal_event}.{contig}.neo.fa",
        wt_fasta="results/microphaser/fasta/{group}/tumor.{tumor_event}.{normal_event}.{contig}.normal.fa",
        tsv="results/microphaser/info/{group}/tumor.{tumor_event}.{normal_event}.{contig}.tsv",
    log:
        "logs/microphaser_tumor/{group}/{tumor_event}.{normal_event}.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser somatic {input.bam} --variants {input.vcf} --ref {input.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.window_length} "
        "< {input.track} > {output.mt_fasta} 2> {log}"


rule microphaser_normal:
    input:
        vcf="results/strelka/normal/{group}.{normal_event}.variants.reheader.norm.bcf",
        bam=get_normal_bam(),
        bai=get_normal_bam(ext=".bam.bai"),
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        wt_fasta=("results/microphaser/fasta/{group}/normal.{normal_event}.{contig}.fa"),
        wt_tsv=("results/microphaser/info/{group}/normal.{normal_event}.{contig}.tsv"),
    log:
        "logs/microphaser_germline/{group}/{normal_event}-{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser normal {input.bam} --variants {input.vcf} --ref {input.ref} -t {output.wt_tsv} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta} 2> {log}"


rule concat_normal_proteome:
    input:
        expand(
            "results/microphaser/fasta/{{group}}/normal.{{normal_event}}.{contig}.fa",
            contig=contigs,
        ),
    output:
        "results/microphaser/fasta/{group}.{normal_event}.normal_proteome.fa",
    log:
        "logs/microphaser/concat_normal_proteome/{group}.{normal_event}.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule build_normal_proteome_db:
    input:
        "results/microphaser/fasta/{group}.{normal_event}.normal_proteome.fa",
    output:
        bin="results/microphaser/bin/{group}.{normal_event}.{mhc}.normal_proteome.bin",
        fasta="results/microphaser/fasta/{group}.{normal_event}.{mhc}.normal_proteome.peptides.fasta",
    log:
        "logs/microphaser/build_normal_proteome_db/{group}.{normal_event}-{mhc}.log",
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
        tsv="results/microphaser/info/{group}/tumor.{tumor_event}.{contig}.tsv",
        proteome=get_proteome(),
    output:
        mt_fasta=(
            "results/microphaser/fasta/filtered/{group}/{mhc}.{tumor_event}.{contig}.neo.fa"
        ),
        wt_fasta=(
            "results/microphaser/fasta/filtered/{group}/{mhc}.{tumor_event}.{contig}.normal.fa"
        ),
        tsv="results/microphaser/info/filtered/{group}/{mhc}.{tumor_event}.{contig}.tsv",
        removed="results/microphaser/info/removed/{group}/{mhc}.{tumor_event}.{contig}.removed.tsv",
    log:
        "logs/microphaser_filter/{group}/{mhc}.{tumor_event}.{contig}.log",
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
            "results/microphaser/info/filtered/{{group}}/{{mhc}}.{{tumor_event}}.{contig}.tsv",
            contig=contigs,
        ),
    output:
        "results/microphaser/info/filtered/{group}.{mhc}.{tumor_event}.tsv",
    log:
        "logs/concat_tsvs/{group}.{mhc}.{tumor_event}.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output} 2>{log}"
