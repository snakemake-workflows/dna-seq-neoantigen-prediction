rule vg2svg:
    input:
        "{prefix}.vl.json"
    output:
        report("{prefix}.svg", caption="../report/tmb.rst", category="Tumor Mutational Burden")
    log:
        "logs/vega/{prefix}.log"
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input} {output} > {log} 2>&1"