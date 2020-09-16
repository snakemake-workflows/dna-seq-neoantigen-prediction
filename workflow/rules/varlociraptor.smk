rule render_scenario:
    input:
        config["calling"]["scenario"]
    output:
        report("results/scenarios/{pair}.yaml", caption="../report/scenario.rst", category="Variant calling scenarios")
    params:
        samples=samples
    log:
        "logs/scenarious/{pair}.log"
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"

rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        candidates="results/candidate-calls/{pair}.{caller}.bcf",
        bcf_index = "results/candidate-calls/{pair}.{caller}.bcf.csi",
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai"
    output:
        "results/observations/{pair}/{sample}.{caller}.bcf"
    log:
        "logs/varlociraptor/preprocess/{pair}/{sample}.{caller}.log"
    conda:
        "../envs/varlociraptor.yaml"
    threads: 8
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} --threads {threads} 2> {log}"
        

rule sort_observations:
    input:
       "results/observations/{pair}/{sample}.{caller}.bcf"
    output:
        "results/observations/{pair}/{sample}.{caller}.sorted.bcf"
    log:
        "logs/sort-observations/{pair}.{sample}.{caller}.log"
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs=get_pair_observations,
        scenario="results/scenarios/{pair}.yaml"
    output:
        temp("results/calls/{pair}.{caller}.bcf")
    log:
        "logs/varlociraptor/call/{pair}.{caller}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_pair_aliases(w), input.obs)]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

rule bcftools_concat:
    input:
        calls = expand(
            "results/calls/{{pair}}.{caller}.bcf",
            caller=caller
        ),
        indexes = expand(
            "results/calls/{{pair}}.{caller}.bcf.csi",
            caller=caller
        ),
    output:
        "results/calls/{pair}.bcf"
    log:
        "logs/bcftools-concat/{pair}.log"
    params:
        "-a -Ob" # Check this
    wrapper:
        "0.36.0/bio/bcftools/concat"
