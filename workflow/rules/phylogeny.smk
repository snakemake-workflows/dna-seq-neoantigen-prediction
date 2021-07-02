def get_somatic_calls(wildcards):
    return expand(
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf",
        sample=samples[samples.type == "tumor"]["sample"],
    )


rule merge_snvs:
    input:
        calls=get_somatic_calls,
    output:
        "results/strelka/merged_calls.vcf",
    log:
        "results/logs/bcftools/merge.log",
    params:
        "--use-header strelka/sampleheader.txt --force-samples",
    wrapper:
        "0.36.0/bio/bcftools/merge"


rule query:
    input:
        "results/strelka/merged_calls.vcf",
    output:
        "results/strelka/call_matrix.tsv",
    log:
        "results/logs/bcftools/query.log",
    params:
        "-H -f '%CHROM\t%POS[\t%DP]\n'",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools query {params} {input} -o {output}"  ## TODO: wrapper


rule nj_tree:
    input:
        matrix="results/strelka/call_matrix.tsv",
    output:
        pdf="results/plots/phylogeny_njtree.pdf",
    log:
        "results/logs/njtree.log",
    conda:
        "../envs/phylogeny.yaml"
    script:
        "../scripts/phylogeny.R"
