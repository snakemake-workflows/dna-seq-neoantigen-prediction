rule annotate_strelka_variants:
    input:
        calls="results/strelka/{calls}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
    output:
        calls="results/strelka/{calls}.annotated.bcf",
        stats=report(
            "results/strelka/{calls}.annotated.stats.html",
            caption="../report/stats.rst",
            category="QC",
        ),
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=[],  #config["annotations"]["vep"]["plugins"],
        extra="--vcf_info_field ANN --hgvs --symbol --canonical",
    log:
        "logs/vep/{calls}.strelka.annotate.log",
    wrapper:
        "0.59.2/bio/vep/annotate"
