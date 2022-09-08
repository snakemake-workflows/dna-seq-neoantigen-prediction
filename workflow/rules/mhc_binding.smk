rule net_mhc_pan:
    input:
        peptides=expand(
            "results/microphaser/fasta/filtered/contigs/{{group}}.{{tumor_alias}}.merged_tumor_normal.pep_len_{peptide_length}.{{contig}}.{{peptide_type}}.fa",
            peptide_length=config["params"]["net_mhc_pan"]["peptide_len"],
        ),
        alleles=get_alleles_MHCI,
    output:
        "results/net_mhc_pan/contigs/{group}.{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.tsv",
    log:
        "logs/net_mhc_pan/{group}.{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.log",
    conda:
        "../envs/tcsh.yaml"
    params:
        extra=config["params"]["net_mhc_pan"]["extra"],
        netMHC=config["params"]["net_mhc_pan"]["location"],
        length=config["params"]["net_mhc_pan"]["peptide_len"],
        alleles=lambda wc, input: ",".join(
            pd.read_csv(input.alleles[0], header=None)[0]
        ),
    shell:
        "( "
        "if [ -s {input.peptides} ]; "
        "then "
        "  {params.netMHC}/netMHCpan {params.extra} -BA -s -l {params.length} -xls -xlsfile {output} -a {params.alleles} -f {input.peptides} > {log}; "
        "else "
        "  touch {output}; "
        "fi "
        " ) 2> {log}"


rule net_mhc_two_pan:
    input:
        peptides=expand(
            "results/microphaser/fasta/filtered/contigs/{{group}}.{{tumor_alias}}.merged_tumor_normal.pep_len_{peptide_length}.{{contig}}.{{peptide_type}}.fa",
            peptide_length=config["params"]["net_mhc_two_pan"]["peptide_len"],
        ),
        alleles=get_alleles_MHCII,
    output:
        "results/net_mhc_two_pan/contigs/{group}.{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.tsv",
    log:
        "logs/net_mhc_two_pan/{group}.{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.log",
    conda:
        "../envs/tcsh.yaml"
    params:
        extra=config["params"]["net_mhc_two_pan"]["extra"],
        netMHC=config["params"]["net_mhc_two_pan"]["location"],
        length=config["params"]["net_mhc_two_pan"]["peptide_len"],
        alleles=lambda wc, input: ",".join(
            pd.read_csv(input.alleles[0], header=None)[0]
        ),
    shell:
        "( "
        "if [ -s {input.peptides} ]; "
        "then "
        "  {params.netMHC}/netMHCIIpan {params.extra} -BA -s -length {params.length} -xls -xlsfile {output} -a {params.alleles} -f {input.peptides} > {log}; "
        "else "
        "  touch {output}; "
        "fi "
        " ) 2> {log}"


rule tidy_mhc_out:
    input:
        expand(
            "results/{{mhc}}/contigs/{{group}}.{{tumor_alias}}.merged_tumor_normal.{contig}.{{peptide_type}}.tsv",
            contig=contigs,
        ),
    output:
        joined_mhc_out="results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.{peptide_type}.tsv",
    log:
        "logs/parse_mhc_out/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.{peptide_type}.log",
    script:
        "../scripts/tidy_mhc_output.py"


rule merge_neoantigen_info:
    input:
        info="results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.tsv",
        neo="results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.neo.tsv",
        normal="results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.normal.tsv",
    output:
        report(
            "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.DNA.tsv",
            caption="../report/neoantigens.dna.rst",
            category="Neoantigens",
        ),
    log:
        "logs/mhc_csv_table/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.log",
    script:
        "../scripts/merge_neoantigen_info.py"


rule add_rna_info:
    input:
        counts="results/kallisto/{group}.{tumor_alias}",
        table="results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.DNA.tsv",
    output:
        report(
            "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.RNA.tsv",
            caption="../report/neoantigens.rna.rst",
            category="Neoantigens",
        ),
    params:
        abundance=lambda wc, input: "{}/abundance.tsv".format(input.counts),
    log:
        "logs/add_rna/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.log",
    script:
        "../scripts/add_rna_info.py"