rule vg2svg:
    input:
        "{prefix}.vl.json"
    output:
        report("{prefix}.svg", caption="../report/tmb.rst", category="Tumor Mutational Burden")
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input} {output}"