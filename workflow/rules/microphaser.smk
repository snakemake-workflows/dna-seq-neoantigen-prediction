rule microphaser_somatic:
    input:
        vcf="strelka/merged/{sample}/all_variants.prepy.bcf",
        bam="bwa/{sample}.rmdup.bam",
        bai="bwa/{sample}.rmdup.bam.bai",
        track=config["reference"]["gtfs"] + "/{chrom}.gtf",
    output:
        mt_fasta="microphaser/fasta/{sample}/{sample}.{chrom}.mt.fa",
        wt_fasta="microphaser/fasta/{sample}/{sample}.{chrom}.wt.fa",
        tsv=temp("microphaser/info/{sample}/{sample}.{chrom}.tsv")
    params:
        ref=config["reference"]["genome"],
        window_length=config["params"]["microphaser"]["window_len"]
    shell:
        "../microphaser/target/release/microphaser somatic {input.bam} --variants {input.vcf} --ref {params.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.window_length} "
        "< {input.track} > {output.mt_fasta}"

rule microphaser_germline:
    input:
        vcf="strelka/germline/{normal}/results/variants/variants.reheader.prepy.bcf",
        bam="bwa/{normal}.rmdup.bam",
        bai="bwa/{normal}.rmdup.bam.bai",
        track=config["reference"]["gtfs"] + "/{chrom}.gtf",
    output:
        wt_fasta="microphaser/fasta/germline/{normal}/{normal}.germline.{chrom}.fa"
    params:
        ref=config["reference"]["genome"],
        window_length=config["params"]["microphaser"]["window_len"]
    shell:
        "../microphaser/target/release/microphaser normal {input.bam} --variants {input.vcf} --ref {params.ref} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta}"

rule concat_proteome:
    input:
        expand("microphaser/fasta/germline/{{normal}}/{{normal}}.germline.{chrom}.fa", chrom = contigs)
    output:
        "microphaser/fasta/germline/{normal}/reference_proteome.fa"
    shell:
        "cat {input} > {output}"

rule build_germline_proteome:
    input:
        "microphaser/fasta/germline/{normal}/reference_proteome.fa"
    output:
        "microphaser/fasta/germline/{normal}/reference_proteome.bin"
    shell:
        "../microphaser/target/release/microphaser build_reference -r {input} -o {output}"

rule microphaser_filter:
    input:
        tsv="microphaser/info/{sample}/{sample}.{chrom}.tsv",
        proteome=get_proteome
    output:
        mt_fasta="microphaser/fasta/{sample}/filtered/{sample}.{chrom}.mt.fa",
        wt_fasta="microphaser/fasta/{sample}/filtered/{sample}.{chrom}.wt.fa",
        tsv="microphaser/info/{sample}/filtered/{sample}.{chrom}.tsv"
    shell:
        "../microphaser/target/release/microphaser filter -r {input.proteome} -t {input.tsv} -o {output.tsv} -n {output.wt_fasta} > {output.mt_fasta}"

rule concat_tsvs:
    input:
        expand("microphaser/info/{{sample}}/filtered/{{sample}}.{chrom}.tsv", chrom = contigs)
    output:
       "microphaser/info/{sample}/filtered/{sample}.tsv"
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output}"

