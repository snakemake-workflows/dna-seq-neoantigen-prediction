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
        "v1.12.0/bio/bcftools/norm"


rule add_somatic_flag:
    input:
        bcf="results/final-calls/{group}.{set}.norm.bcf",
        header_line="resources/somatic_flag_header_line.txt",
        flag_bed="resources/genome.somatic_flag.bed.gz",
        flag_bed_idx="resources/genome.somatic_flag.bed.gz.tbi",
    output:
        "results/final-calls/{group}.{set}.somatic_flag.norm.bcf",
    log:
        "logs/bcftools_annotate/{group}.{set}.somatic_flag.norm.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bcftools annotate "
        "  --annotations {input.flag_bed} "
        "  --mark-sites +SOMATIC "
        "  --columns CHROM,FROM,TO "
        "  --header-lines {input.header_line} "
        "  -O b -o {output} "
        "  {input.bcf} "
        ") 2> {log}"


rule merge_tumor_normal:
    input:
        calls=expand(
            "results/final-calls/{{group}}.{sets}.norm.bcf",
            sets=[
                config["params"]["microphaser"]["events"]["normal"],
                config["params"]["microphaser"]["events"]["tumor"] + ".somatic_flag",
            ],
        ),
        index=expand(
            "results/final-calls/{{group}}.{sets}.norm.bcf.csi",
            sets=[
                config["params"]["microphaser"]["events"]["normal"],
                config["params"]["microphaser"]["events"]["tumor"] + ".somatic_flag",
            ],
        ),
    output:
        "results/final-calls/{group}.merged_tumor_normal.norm.bcf",
    log:
        "logs/bcftools/concat-tumor-normal/{group}.merged_tumor_normal.log",
    params:
        extra="-O b -a",
    wrapper:
        "v1.12.0/bio/bcftools/concat"


rule microphaser_tumor:
    input:
        bcf="results/final-calls/{group}.merged_tumor_normal.norm.annotated.bcf",
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        track="resources/annotation/{contig}.gtf",
        ref="resources/genome.fasta",
    output:
        mt_fasta="results/microphaser/fasta/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.neo.fa",
        wt_fasta="results/microphaser/fasta/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.normal.fa",
        tsv="results/microphaser/info/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.tsv",
    log:
        "logs/microphaser_tumor/{group}/{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=lambda wc: int(wc.peptide_length) * 3,
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
        wt_fasta=(
            "results/microphaser/fasta/contigs/{group}.{normal_alias}.{normal_set}.pep_len_{peptide_length}.{contig}.fa"
        ),
        wt_tsv=(
            "results/microphaser/info/contigs/{group}.{normal_alias}.{normal_set}.pep_len_{peptide_length}.{contig}.tsv"
        ),
    log:
        "logs/microphaser_normal/contigs/{group}/{normal_alias}.{normal_set}.pep_len_{peptide_length}.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    params:
        window_length=lambda wc: int(wc.peptide_length) * 3,
    shell:
        "microphaser normal {input.bam} --variants {input.bcf} --ref {input.ref} -t {output.wt_tsv} -w {params.window_length} "
        "< {input.track} > {output.wt_fasta} 2> {log}"


rule concat_normal_proteome:
    input:
        expand(
            "results/microphaser/fasta/contigs/{{group}}.normal.{{normal_set}}.pep_len_{{peptide_length}}.{contig}.fa",
            contig=contigs,
        ),
    output:
        "results/microphaser/fasta/{group}.{normal_set}.normal_proteome.pep_len_{peptide_length}.fa",
    log:
        "logs/microphaser_concat_normal_proteome/{group}.{normal_set}.pep_len_{peptide_length}.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule build_normal_proteome_db:
    input:
        "results/microphaser/fasta/{group}.{normal_set}.normal_proteome.pep_len_{peptide_length}.fa",
    output:
        bin="results/microphaser/bin/{group}.{normal_set}.normal_proteome.pep_len_{peptide_length}.bin",
        fasta="results/microphaser/fasta/{group}.{normal_set}.normal_proteome.pep_len_{peptide_length}.peptides.fasta",
    log:
        "logs/microphaser_build_normal_proteome_db/{group}.{normal_set}.pep_len_{peptide_length}.log",
    conda:
        "../envs/microphaser.yaml"
    shell:
        "( microphaser build_reference -r {input} -o {output.bin} -l {wildcards.peptide_length} > {output.fasta} ) 2> {log}"


rule microphaser_filter:
    input:
        tsv="results/microphaser/info/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.tsv",
        proteome=expand(
            "results/microphaser/bin/{{group}}.{normal_set}.normal_proteome.pep_len_{{peptide_length}}.bin",
            normal_set=config["params"]["microphaser"]["events"]["normal"],
        ),
    output:
        mt_fasta=(
            "results/microphaser/fasta/filtered/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.neo.fa"
        ),
        wt_fasta=(
            "results/microphaser/fasta/filtered/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.normal.fa"
        ),
        tsv="results/microphaser/info/filtered/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.tsv",
        removed="results/microphaser/info/removed/contigs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.removed.tsv",
    log:
        "logs/microphaser_filter/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.{contig}.log",
    conda:
        "../envs/microphaser.yaml"
    shell:
        "microphaser filter -r {input.proteome} -t {input.tsv} -o {output.tsv} -n {output.wt_fasta} -s {output.removed} -l {wildcards.peptide_length} > {output.mt_fasta} 2>{log}"


rule concat_tsvs:
    input:
        expand(
            "results/microphaser/info/filtered/contigs/{{group}}.{{tumor_alias}}.merged_tumor_normal.pep_len_{{peptide_length}}.{contig}.tsv",
            contig=contigs,
        ),
    output:
        "results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.tsv",
    log:
        "logs/microphaser_concat_tsvs/{group}.{tumor_alias}.merged_tumor_normal.pep_len_{peptide_length}.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d '\t' {input} | xsv fmt -t '\t' -d ',' > {output} 2>{log}"
