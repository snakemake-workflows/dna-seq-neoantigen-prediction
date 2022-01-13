rule filter_by_annotation:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.{filter}.filtered_ann.bcf",
    log:
        "logs/filter-calls/annotation/{prefix}.{filter}.log",
    params:
        filter=lambda w: config["calling"]["filter"][w.filter],
    conda:
        "../envs/vembrane.yaml"
    shell:
        'vembrane filter --output-fmt bcf --output {output} "{params.filter}" {input} &> {log}'


rule filter_odds:
    input:
        get_annotated_bcf,
    output:
        "results/calls/{cancer_sample}.{event}.{scatteritem}.filtered_odds.bcf",
    params:
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event][
            "varlociraptor"
        ],
    log:
        "logs/filter-calls/posterior_odds/{cancer_sample}.{scatteritem}.{event}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule gather_calls:
    input:
        calls=gather.calling(
            "results/calls/{{cancer_sample}}.{{event}}.{scatteritem}.filtered_odds.bcf"
        ),
        idx=gather.calling(
            "results/calls/{{cancer_sample}}.{{event}}.{scatteritem}.filtered_odds.bcf.csi"
        ),
    output:
        "results/calls/{cancer_sample}.{event}.filtered_odds.bcf",
    log:
        "logs/gather-calls/{cancer_sample}.{event}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.67.0/bio/bcftools/concat"


rule control_fdr:
    input:
        "results/calls/{cancer_sample}.{event}.filtered_odds.bcf",
    output:
        "results/calls/{cancer_sample}.{vartype}.{event}.fdr-controlled.bcf",
    log:
        "logs/control-fdr/{cancer_sample}.{vartype}.{event}.log",
    params:
        query=get_fdr_control_params,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.query[events]} --fdr {params.query[threshold]} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_input(".bcf"),
        idx=get_merge_input(".bcf.csi"),
    output:
        "results/merged-calls/{cancer_sample}.{event}.fdr-controlled.bcf",
    log:
        "logs/merge-calls/{cancer_sample}.{event}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.59.2/bio/bcftools/concat"


rule change_samplenames:
    input:
        call="results/merged-calls/{cancer_sample}.{event}.fdr-controlled.bcf",
    output:
        temp("results/merged-calls/{cancer_sample}.{event}.renaming.txt"),
    log:
        "logs/change-samplenames/{cancer_sample}.{event}.log",
    params:
        prefix=lambda w, input: os.path.basename(input["call"]).split(".")[0],
    shell:
        "echo -e 'normal {params.prefix}_N\ntumor {params.prefix}_T' > {output}"


rule reheader_varlociraptor:
    input:
        vcf="results/merged-calls/{cancer_sample}.{event}.fdr-controlled.bcf",
        samples="results/merged-calls/{cancer_sample}.{event}.renaming.txt",
    output:
        "results/merged-calls/{cancer_sample}.{event}.reheader.bcf",
    log:
        "logs/reheader-calls/{cancer_sample}.{event}.log",
    params:
        extra="",
        view_extra="-O b",
    wrapper:
        "0.60.0/bio/bcftools/reheader"
