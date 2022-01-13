rule freebayes:
    input:
        ref="resources/genome.fasta",
        # you can have a list of samples here
        samples=get_paired_bams,
    output:
        "results/candidate-calls/{cancer_sample}.freebayes.bcf",
    log:
        "logs/{cancer_sample}.log",
    params:
        extra=config["params"].get("freebayes", ""),
        chunksize=100000,
    threads: 60
    wrapper:
        "0.65.0/bio/freebayes"


rule scatter_candidates:
    input:
        "results/candidate-calls/{cancer_sample}.{caller}.bcf",
    output:
        scatter.calling("results/candidate-calls/{{cancer_sample}}.{{caller}}.{scatteritem}.bcf"),
    log:
        "logs/scatter-candidates/{cancer_sample}.{caller}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"
