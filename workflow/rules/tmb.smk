if config["tmb"]["activate"]:
    rule estimate_tmb:
        input:
            "results/calls/{pair}.annotated.bcf"
        output:
            "results/plots/tmb/{pair}.{plotmode}.vl.json",
        conda:
            "../envs/varlociraptor.yaml"
        params:
            **config["tmb"]
        shell:
            "varlociraptor estimate tmb "
            " --plot-mode {wildcards.plotmode} "
            "--coding-genome-size {params.coding_genome_size} "
            "--somatic-tumor-events {params.somatic_events} "
            "--tumor-sample {params.tumor_sample} "
            "< {input} > {output}"
