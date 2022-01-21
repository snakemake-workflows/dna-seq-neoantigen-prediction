rule microphaser_somatic:
    input:
        vcf="results/strelka/merged/{cancer_sample}/all_variants.norm.annotated.bcf",
        bam="results/recal/{cancer_sample}.sorted.bam",
        bai="results/recal/{cancer_sample}.sorted.bam.bai",
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        mt_fasta="results/microphaser/fasta/{cancer_sample}/{cancer_sample}.{contig}.neo.fa",
        wt_fasta="results/microphaser/fasta/{cancer_sample}/{cancer_sample}.{contig}.normal.fa",
        tsv="results/microphaser/info/{cancer_sample}/{cancer_sample}.{contig}.tsv",
    log:
        "logs/microphaser/somatic/{cancer_sample}-{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser somatic {input.bam} --variants {input.vcf} --ref {input.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.window_length} "
        "< {input.track} > {output.mt_fasta} 2> {log}"


rule microphaser_germline:
    input:
        vcf="results/strelka/germline/{normal_sample}/results/variants/variants.reheader.norm.bcf",
        bam="results/recal/{normal_sample}.sorted.bam",
        bai="results/recal/{normal_sample}.sorted.bam.bai",
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        wt_fasta=(
            "results/microphaser/fasta/germline/{normal_sample}/{normal_sample}.germline.{contig}.fa"
        ),
        wt_tsv=(
            "results/microphaser/info/germline/{normal_sample}/{normal_sample}.germline.{contig}.tsv"
        ),
    log:
        "logs/microphaser/germline/{normal_sample}-{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=config["params"]["microphaser"]["window_len"],
    shell:
        "microphaser normal {input.bam} --variants {input.vcf} --ref {input.ref} -t {output.wt_tsv} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta} 2> {log}"


rule concat_proteome:
    input:
        expand(
            "results/microphaser/fasta/germline/{{normal_sample}}/{{normal_sample}}.germline.{contig}.fa",
            contig=contigs,
        ),
    output:
        "results/microphaser/fasta/germline/{normal_sample}/reference_proteome.fa",
    log:
        "logs/microphaser/concat-ref-proteome/{normal_sample}.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule build_germline_proteome:
    input:
        "results/microphaser/fasta/germline/{normal_sample}/reference_proteome.fa",
    output:
        bin="results/microphaser/fasta/germline/{normal_sample}/{mhc}/reference_proteome.bin",
        fasta="results/microphaser/fasta/germline/{normal_sample}/{mhc}/reference_proteome.peptides.fasta",
    log:
        "logs/microphaser/build-ref-proteome-db/{normal_sample}-{mhc}.log",
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
        tsv="results/microphaser/info/{cancer_sample}/{cancer_sample}.{contig}.tsv",
        proteome=get_proteome,
    output:
        mt_fasta=(
            "results/microphaser/fasta/{cancer_sample}/filtered/{mhc}/{cancer_sample}.{contig}.neo.fa"
        ),
        wt_fasta=(
            "results/microphaser/fasta/{cancer_sample}/filtered/{mhc}/{cancer_sample}.{contig}.normal.fa"
        ),
        tsv="results/microphaser/info/{cancer_sample}/filtered/{mhc}/{cancer_sample}.{contig}.tsv",
        removed="results/microphaser/info/{cancer_sample}/removed/{mhc}/{cancer_sample}.{contig}.removed.tsv",
    log:
        "logs/microphaser/filter/{cancer_sample}-{mhc}-{contig}.log",
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
            "results/microphaser/info/{{cancer_sample}}/filtered/{{mhc}}/{{cancer_sample}}.{contig}.tsv",
            contig=contigs,
        ),
    output:
        "results/microphaser/info/{cancer_sample}/filtered/{mhc}/{cancer_sample}.tsv",
    log:
        "logs/concat-tsv/{cancer_sample}-{mhc}.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output} 2>{log}"
