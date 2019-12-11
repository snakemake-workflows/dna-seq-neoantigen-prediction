rule render_scenario:
    input:
        config["calling"]["scenario"]
    output:
        report("scenarios/{pair}.yaml", caption="../report/scenario.rst", category="Variant calling scenarios")
    params:
        samples = samples
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"

rule varlociraptor_preprocess:
    input:
        ref=config["reference"]["genome"],
        candidates="candidate-calls/{pair}.{caller}.bcf",
        bcf_index = "candidate-calls/{pair}.{caller}.bcf.csi",
        bam="bwa/{sample}.rmdup.bam",
        bai="bwa/{sample}.rmdup.bam.bai"
    output:
        "observations/{pair}/{sample}.{caller}.{contig}.bcf"
    log:
        "logs/varlociraptor/preprocess/{pair}/{sample}.{caller}.{contig}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "bcftools view -Ou {input.candidates} {wildcards.contig} | "
        "varlociraptor preprocess variants "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs=get_pair_observations,
        scenario="scenarios/{pair}.yaml"
    output:
        temp("calls/{pair}.{caller}.{contig}.bcf")
    log:
        "logs/varlociraptor/call/{pair}.{caller}.{contig}.log"
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
            "calls/{{pair}}.{caller}.{contig}.bcf",
            caller=caller,
            contig=contigs
        ),
        indexes = expand(
            "calls/{{pair}}.{caller}.{contig}.bcf.csi",
            caller=caller,
            contig=contigs
        ),
    output:
        "calls/{pair}.vcf"
    params:
        "-a" # Check this
    wrapper:
        "0.36.0/bio/bcftools/concat"
