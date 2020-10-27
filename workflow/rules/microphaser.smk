rule microphaser_somatic:
    input:
        vcf="results/strelka/merged/{sample}/all_variants.norm.annotated.bcf",
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai",
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta"
    output:
        mt_fasta="results/microphaser/fasta/{sample}/{sample}.{contig}.mt.fa",
        wt_fasta="results/microphaser/fasta/{sample}/{sample}.{contig}.wt.fa",
        tsv="results/microphaser/info/{sample}/{sample}.{contig}.tsv"
    log:
        "logs/microphaser/somatic/{sample}-{contig}.log"
    params:
        window_length=config["params"]["microphaser"]["window_len"]
    shell:
        "../microphaser/target/release/microphaser somatic {input.bam} --variants {input.vcf} --ref {input.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.window_length} "
        "< {input.track} > {output.mt_fasta} 2> {log}"

rule microphaser_germline:
    input:
        vcf="results/strelka/germline/{normal}/results/variants/variants.reheader.prepy.bcf",
        bam="results/recal/{normal}.sorted.bam",
        bai="results/recal/{normal}.sorted.bam.bai",
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta"
    output:
        wt_fasta="results/microphaser/fasta/germline/{normal}/{normal}.germline.{contig}.fa",
        wt_tsv="results/microphaser/info/germline/{normal}/{normal}.germline.{contig}.tsv"
    log:
        "logs/microphaser/germline/{normal}-{contig}.log"
    params:
        window_length=config["params"]["microphaser"]["window_len"]
    shell:
        "../microphaser/target/release/microphaser normal {input.bam} --variants {input.vcf} --ref {input.ref} -t {output.wt_tsv} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta} 2> {log}"

rule concat_proteome:
    input:
        expand("results/microphaser/fasta/germline/{{normal}}/{{normal}}.germline.{contig}.fa", contig = contigs)
    output:
        "results/microphaser/fasta/germline/{normal}/reference_proteome.fa"
    log:
        "logs/microphaser/concat-ref-proteome/{normal}.log"
    shell:
        "cat {input} > {output} 2> {log}"

rule build_germline_proteome:
    input:
        "results/microphaser/fasta/germline/{normal}/reference_proteome.fa"
    output:
        bin="results/microphaser/fasta/germline/{normal}/{mhc}/reference_proteome.bin",
        fasta="results/microphaser/fasta/germline/{normal}/{mhc}/reference_proteome.peptides.fasta"
    log:
        "logs/microphaser/build-ref-proteome-db/{normal}-{mhc}.log"
    params:
        length=lambda wildcards: config["params"]["microphaser"]["peptide_len"][wildcards.mhc]
    shell:
        "../microphaser/target/release/microphaser build_reference -r {input} -o {output.bin} -l {params.length} --peptides {output.fasta} > {log} 2>&1"

rule microphaser_filter:
    input:
        tsv="results/microphaser/info/{sample}/{sample}.{contig}.tsv",
        proteome=get_proteome
    output:
        mt_fasta="results/microphaser/fasta/{sample}/filtered/{mhc}/{sample}.{contig}.mt.fa",
        wt_fasta="results/microphaser/fasta/{sample}/filtered/{mhc}/{sample}.{contig}.wt.fa",
        tsv="results/microphaser/info/{sample}/filtered/{mhc}/{sample}.{contig}.tsv",
        removed="results/microphaser/info/{sample}/removed/{mhc}/{sample}.{contig}.removed.tsv"
    log:
        "logs/microphaser/filter/{sample}-{mhc}-{contig}.log"
    params:
        length=lambda wildcards: config["params"]["microphaser"]["peptide_len"][wildcards.mhc]
    shell:
        "../microphaser/target/release/microphaser filter -r {input.proteome} -t {input.tsv} -o {output.tsv} -n {output.wt_fasta} -s {output.removed} -l {params.length} > {output.mt_fasta} 2>{log}"


rule concat_tsvs:
    input:
        expand("results/microphaser/info/{{sample}}/filtered/{{mhc}}/{{sample}}.{contig}.tsv", contig = contigs)
    output:
       "results/microphaser/info/{sample}/filtered/{mhc}/{sample}.tsv"
    log:
        "logs/concat-tsv/{sample}-{mhc}.log"
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output} 2>{log}"
