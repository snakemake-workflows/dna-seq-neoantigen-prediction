# rule filter_by_annotation:
#     input:
#         get_annotated_bcf
#     output:
#         "results/calls/{group}.{filter}.filtered_ann.bcf"
#     log:
#         "logs/filter-calls/annotation/{group}.{filter}.log"
#     params:
#         filter=lambda w: config["calling"]["filter"][w.filter]
#     conda:
#         "../envs/vembrane.yaml"
#     shell:
#         "vembrane {params.filter:q} {input} --output-fmt bcf --output {output} &> {log}"


rule filter_odds:
    input:
        get_annotated_bcf
    output:
        "results/calls/{pair}.{event}.filtered_odds.bcf"
    params:
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    log:
        "logs/filter-calls/posterior_odds/{pair}.{event}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"



rule control_fdr:
    input:
        "results/calls/{pair}.{event}.filtered_odds.bcf"
    output:
        "results/calls/{pair}.{vartype}.{event}.fdr-controlled.bcf"
    params:
        threshold=config["calling"]["fdr-control"]["threshold"],
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    log:
        "logs/control-fdr/{pair}-{vartype}-{event}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.threshold} > {output}"


rule merge_calls:
    input:
        calls=get_merge_input(".bcf"),
        idx=get_merge_input(".bcf.csi")
    output:
        "results/merged-calls/{pair}.{event}.fdr-controlled.bcf"
    log:
        "logs/merge-calls/{pair}-{event}.log"
    params:
        "-a -Ob"
    wrapper:
        "0.37.1/bio/bcftools/concat"
