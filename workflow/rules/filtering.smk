rule control_fdr:
    input:
        "calls/{group}.vcf"
    output:
        "calls/{group}.{event}.{vartype}.fdr-controlled.bcf"
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
        calls=expand("calls/{{group}}.{{event}}.{vartype}.fdr-controlled.bcf", vartype=["SNV", "DEL", "INS"]),
        indexes=expand("calls/{{group}}.{{event}}.{vartype}.fdr-controlled.bcf.csi", vartype=["SNV", "DEL", "INS"])
    output:
        "calls/{group}.{event}.fdr-controlled.vcf"
    params:
        "-a" # Check this
    wrapper:
        "0.35.2/bio/bcftools/concat"
