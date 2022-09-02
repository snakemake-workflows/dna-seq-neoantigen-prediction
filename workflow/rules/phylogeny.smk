def get_somatic_calls(wildcards):
    return expand(
        "results/final-calls/somatic/{sample}/results/variants/somatic.complete.tumor.bcf",
        sample=samples[samples.alias == "tumor"]["sample_name"],
    )


rule merge_snvs:
    input:
        calls=get_somatic_calls,
    output:
        "results/final-calls/merged_calls.vcf",
    log:
        "results/logs/bcftools/merge.log",
    params:
        "--use-header final-calls/sampleheader.txt --force-samples",
    wrapper:
        "v1.12.0/bio/bcftools/merge"


rule query:
    input:
        "results/final-calls/merged_calls.vcf",
    output:
        "results/final-calls/call_matrix.tsv",
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
        matrix="results/final-calls/call_matrix.tsv",
    output:
        pdf="results/plots/phylogeny_njtree.pdf",
    log:
        "results/logs/njtree.log",
    conda:
        "../envs/phylogeny.yaml"
    script:
        "../scripts/phylogeny.R"
