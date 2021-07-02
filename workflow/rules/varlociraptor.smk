rule render_scenario:
    input:
        config["calling"]["scenario"],
    output:
        report(
            "results/scenarios/{pair}.yaml",
            caption="../report/scenario.rst",
            category="Variant calling scenarios",
        ),
    params:
        samples=samples,
    log:
        "logs/scenarious/{pair}.log",
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates="results/candidate-calls/{pair}.{caller}.{scatteritem}.bcf",
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
    output:
        "results/observations/{pair}/{sample}.{caller}.{scatteritem}.bcf",
    params:
        omit_isize="",
    log:
        "logs/varlociraptor/preprocess/{pair}/{sample}.{caller}.{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.omit_isize} --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs=get_pair_observations,
        scenario="results/scenarios/{pair}.yaml",
    output:
        temp("results/calls/{pair}.{caller}.{scatteritem}.bcf"),
    log:
        "logs/varlociraptor/call/{pair}.{caller}.{scatteritem}.log",
    params:
        obs=lambda w, input: [
            "{}={}".format(s, f) for s, f in zip(get_pair_aliases(w), input.obs)
        ],
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
        "benchmarks/varlociraptor/call/{pair}.{caller}.{scatteritem}.tsv"
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"


rule sort_calls:
    input:
        "results/calls/{pair}.{caller}.{scatteritem}.bcf",
    output:
        temp("results/calls/{pair}.{caller}.{scatteritem}.sorted.bcf"),
    log:
        "logs/bcf-sort/{pair}.{caller}.{scatteritem}.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000,
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_concat:
    input:
        calls=get_scattered_calls(),
        indexes=get_scattered_calls(ext=".bcf.csi"),
    output:
        "results/calls/{pair}.{scatteritem}.bcf",
    log:
        "logs/concat-calls/{pair}.{scatteritem}.log",
    params:
        "-a -Ob",  # TODO Check this
    wrapper:
        "0.59.2/bio/bcftools/concat"
