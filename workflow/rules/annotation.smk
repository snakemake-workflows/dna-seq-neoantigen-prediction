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
