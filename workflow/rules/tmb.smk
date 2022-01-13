if config["tmb"]["activate"]:

    rule estimate_tmb:
        input:
            "results/merged-calls/{cancer_sample}.somatic.fdr-controlled.bcf",
        output:
            "results/plots/tmb/{cancer_sample}.{plotmode}.vl.json",
        conda:
            "../envs/varlociraptor.yaml"
        log:
            "logs/tmb/{cancer_sample}-{plotmode}.log",
        params:
            **config["tmb"],
        shell:
            "varlociraptor estimate tmb "
            " --plot-mode {wildcards.plotmode} "
            "--coding-genome-size {params.coding_genome_size} "
            "--somatic-tumor-events {params.somatic_events} "
            "--tumor-sample {params.tumor_sample} "
            "< {input} > {output} 2> {log}"
