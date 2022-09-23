rule prepare_neoprint:
    input:
        neopeptides="results/neo_fox/annotated/{group}.{tumor_alias}.annotated_neoantigens.tsv",
    output:
        mhc_one="results/tables/neoprint/{group}.{tumor_alias}.annotated_neopeptides.I.sorted.tsv",
        mhc_two="results/tables/neoprint/{group}.{tumor_alias}.annotated_neopeptides.II.sorted.tsv",
    log:
        "logs/prepare_neoprint/{group}.{tumor_alias}.log",
    params:
        purity = lambda wc: samples.loc[(samples["group"] == wc.group) & (samples["alias"] == wc.tumor_alias), "purity"].squeeze()
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/prepare_neoprint.py"


rule render_datavzrd_neoprint_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/neo_fox-neoantigens-template.datavzrd.yaml"
        ),
        neopeptides="results/tables/neoprint/{group}.{tumor_alias}.annotated_neopeptides.{mhc}.sorted.tsv",
    output:
        "resources/datavzrd/{group}.{tumor_alias}.datavzrd_neoprint.{mhc}.yaml",
    log:
        "logs/datavzrd_render_neoprint/{group}.{tumor_alias}.{mhc}.log",
    template_engine:
        "yte"


rule datavzrd_neoprint:
    input:
        neopeptides="results/tables/neoprint/{group}.{tumor_alias}.annotated_neopeptides.{mhc}.sorted.tsv",
        config="resources/datavzrd/{group}.{tumor_alias}.datavzrd_neoprint.{mhc}.yaml",
    output:
        report(
            directory("results/datavzrd/neoprint/{group}.{tumor_alias}.{mhc}"),
            htmlindex="index.html",
            caption="../report/neopeptides.rst",
            category="Neopeptides",
            labels=lambda wc: {"group": wc.group, "sample_type": wc.tumor_alias, "MHC": wc.mhc},
        ),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_neoprint/{group}.{tumor_alias}.{mhc}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
