rule filter_by_annotation:
    input:
        get_annotated_bcf
    output:
        "results/calls/{group}.{filter}.filtered.bcf"
    params:
        filter=lambda w: config["calling"]["filter"][w.filter]
    conda:
        "../envs/snpsift.yaml"
    shell:
        "bcftools view {input} | SnpSift filter \"{params.filter}\" | bcftools view -Ob > {output}"


rule control_fdr:
    input:
        "results/calls/{pair}.vcf"
    output:
        "results/calls/{pair}.{event}.{vartype}.fdr-controlled.bcf"
    params:
        threshold=config["calling"]["fdr-control"]["threshold"],
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.threshold} > {output}"

rule concat_vartypes:
    input:
        calls=expand("results/calls/{{pair}}.{{event}}.{vartype}.fdr-controlled.bcf", vartype=["SNV", "DEL", "INS"]),
        indexes=expand("results/calls/{{pair}}.{{event}}.{vartype}.fdr-controlled.bcf.csi", vartype=["SNV", "DEL", "INS"])
    output:
        "calls/{pair}.{event}.fdr-controlled.vcf"
    params:
        "-a" # Check this
    wrapper:
        "0.35.2/bio/bcftools/concat"
