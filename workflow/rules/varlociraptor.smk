rule render_scenario:
    input:
        config["calling"]["scenario"]
    output:
        report("scenarios/{group}.yaml", caption="../report/scenario.rst", category="Variant calling scenarios")
    params:
        samples = samples
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"

rule varlociraptor_preprocess:
    input:
        ref=config["reference"]["genome"],
        candidates="candidate-calls/{group}.{caller}.bcf",
        bcf_index = "candidate-calls/{group}.{caller}.bcf.csi",
        bam="bwa/{sample}.rmdup.bam",
        bai="bwa/{sample}.rmdup.bam.bai"
    output:
        "observations/{group}/{sample}.{caller}.{contig}.bcf"
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.{contig}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "bcftools view -Ou {input.candidates} {wildcards.contig} | "
        "varlociraptor preprocess variants "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="scenarios/{group}.yaml"
    output:
        temp("calls/{group}.{caller}.{contig}.bcf")
    log:
        "logs/varlociraptor/call/{group}.{caller}.{contig}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

rule bcftools_concat:
    input:
        calls = expand(
            "calls/{{group}}.{caller}.{contig}.bcf",
            caller=caller,
            contig=contigs
        ),
        indexes = expand(
            "calls/{{group}}.{caller}.{contig}.bcf.csi",
            caller=caller,
            contig=contigs
        ),
    output:
        "calls/{group}.vcf"
    params:
        "-a" # Check this
    wrapper:
        "0.36.0/bio/bcftools/concat"
