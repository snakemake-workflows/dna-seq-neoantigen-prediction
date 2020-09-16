rule annotate_variants:
    input:
        calls="results/calls/{group}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins"
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats=report("results/calls/{group}.stats.html", caption="../report/stats.rst", category="QC")
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra="{} --vcf_info_field ANN".format(config["annotations"]["vep"]["params"])
    log:
        "logs/vep/{group}.annotate.log"
    wrapper:
        "0.59.2/bio/vep/annotate"

rule annotate_strelka_variants:
    input:
        calls="results/strelka/merged/{sample}/all_variants.norm.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins"
    output:
        calls="results/strelka/merged/{sample}/all_variants.norm.annotated.bcf",
        stats=report("results/strelka/merged/{sample}/{sample}.annotated.stats.html", caption="../report/stats.rst", category="QC")
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=[],#config["annotations"]["vep"]["plugins"],
        extra="--vcf_info_field ANN --hgvs --symbol --canonical"
    log:
        "logs/vep/{sample}.strelka.annotate.log"
    wrapper:
        "0.59.2/bio/vep/annotate"

