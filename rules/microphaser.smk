rule microphaser_somatic:
    input:
        vcf="strelka/{tumor}-{normal}/results/variants/all_variants.prepy.bcf",
        bam="bwa/{tumor}.rmdup.bam",
        bai="bwa/{tumor}.rmdup.bam.bai",
        track="../gtfs/{chrom}.gtf",
    output:
        mt_fasta="microphaser/fasta/{tumor}-{normal}/{tumor}-{normal}.{chrom}.mt.fa",
        wt_fasta="microphaser/fasta/{tumor}-{normal}/{tumor}-{normal}.{chrom}.wt.fa",
        tsv="microphaser/info/{tumor}-{normal}/{tumor}-{normal}.{chrom}.tsv"
    params:
        ref=config["reference"]["genome"],
        windowlen=config["params"]["microphaser"]["windowlen"]
    shell:
        "../microphaser/target/release/microphaser somatic {input.bam} "
        "--variants {input.vcf} --ref {params.ref} --tsv {output.tsv} -n {output.wt_fasta} -w {params.windowlen} "
        "< {input.track} > {output.mt_fasta}"

rule microphaser_germline:
    input:
        vcf="strelka/{normal}/results/variants/variants.reheader.prepy.bcf",
        bam="bwa/{normal}.rmdup.bam",
        bai="bwa/{normal}.rmdup.bam.bai",
        track="../gtfs/{chrom}.gtf",
    output:
        wt_fasta="microphaser/fasta/{normal}/{normal}.{chrom}.fa"
    params:
        ref=config["reference"]["genome"],
        windowlen=config["params"]["microphaser"]["windowlen"]
    shell:
        "../microphaser/target/release/microphaser normal {input.bam} "
        "--variants {input.vcf} --ref {params.ref} -w {params.windowlen} "
        "< {input.track} > {output.wt_fasta}"

rule concat_proteome:
    input:
        expand("microphaser/fasta/{{normal}}/{{normal}}.{chrom}.fa", chrom = CHROMOSOMES)
    output:
        "microphaser/fasta/{normal}/reference_proteome.fa"
    shell:
        "cat {input} > {output}"

rule build_germline_proteome:
    input:
        "microphaser/fasta/{normal}/reference_proteome.fa"
    output:
        "microphaser/fasta/{normal}/reference_proteome.bin"
    shell:
        "../microphaser/target/release/microphaser build_reference -r {input} -o {output}"

rule microphaser_filter:
    input:
        tsv="microphaser/info/{tumor}-{normal}/{tumor}-{normal}.{chrom}.tsv",
        proteome="microphaser/fasta/{normal}/reference_proteome.bin"
    output:
        mt_fasta="microphaser/fasta/{tumor}-{normal}/filtered/{tumor}-{normal}.{chrom}.mt.fa",
        wt_fasta="microphaser/fasta/{tumor}-{normal}/filtered/{tumor}-{normal}.{chrom}.wt.fa",
        tsv="microphaser/info/{tumor}-{normal}/filtered/{tumor}-{normal}.{chrom}.tsv"
    shell:
        "../microphaser/target/release/microphaser filter "
        "-r {input.proteome} -t {input.tsv} -o {output.tsv} -n {output.wt_fasta} "
        "> {output.mt_fasta}"

rule concat_tsvs:
    input:
        expand("microphaser/info/{{tumor}}-{{normal}}/filtered/{{tumor}}-{{normal}}.{chrom}.tsv", chrom = CHROMOSOMES)
    output:
       "microphaser/info/{tumor}-{normal}/filtered/{tumor}-{normal}.tsv"
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output}"

