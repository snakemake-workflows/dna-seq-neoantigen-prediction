# rule mhcflurry:
#     input:
#         peptides="results/microphaser/fasta/{sample}/filtered/{sample}.{contig}.{peptide_type}.fa",
#         alleles="results/optitype/{sample}/hla_alleles_{sample}.tsv",
#         wt_alleles=get_germline_optitype
#     output:
#         "results/mhcflurry/{sample}/{contig}/output.{peptide_type}.csv"
#     log:
#         "logs/mhcflurry/{sample}-{contig}-{peptide_type}.log"
#     run:
#         if "wt" in input.peptides:
#             alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
#         else:
#             alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
#         cmd = "if [ -s {input.peptides} ]; then mhctools --mhc-predictor mhcflurry --mhc-alleles {alleles} --input-fasta-file {input.peptides} --output-csv {output} > {log}; else touch {output}; fi"
#         shell(cmd)


rule netMHCpan:
    input:
        peptides="results/microphaser/fasta/filtered/{group}/{tumor_alias}.merged_tumor_normal.netMHCpan.{contig}.{peptide_type}.fa",
        alleles=get_alleles_MHCI,
    output:
        "results/netMHCpan/{group}/{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.xls",
    log:
        "logs/netMHCpan/{group}/{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.log",
    conda:
        "../envs/tcsh.yaml"
    params:
        extra=config["affinity"]["netMHCpan"]["params"],
        netMHC=config["affinity"]["netMHCpan"]["location"],
        alleles=lambda wc, input: ",".join(
            pd.read_csv(input.alleles[0], sep="\t").iloc[0]
        ),
    shell:
        "( "
        "if [ -s {input.peptides} ]; "
        "then "
        "  {params.netMHC}/netMHCpan {params.extra} -xlsfile {output} -a {params.alleles} -f {input.peptides} > {log}; "
        "else "
        "  touch {output}; "
        "fi "
        " ) 2> {log}"


rule netMHCIIpan:
    input:
        peptides="results/microphaser/fasta/filtered/{group}/{tumor_alias}.merged_tumor_normal.netMHCIIpan.{contig}.{peptide_type}.fa",
        alleles=get_alleles_MHCII,
    output:
        "results/netMHCIIpan/{group}/{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.xls",
    log:
        "logs/netMHCIIpan/{group}/{tumor_alias}.merged_tumor_normal.{contig}.{peptide_type}.log",
    conda:
        "../envs/tcsh.yaml"
    params:
        extra=config["affinity"]["netMHCIIpan"]["params"],
        netMHC=config["affinity"]["netMHCIIpan"]["location"],
        alleles=lambda wc, input: ",".join(
            pd.read_csv(input.alleles[0], sep="\t")["Allele"].tolist()
        ),
    shell:
        "( "
        "if [ -s {input.peptides} ]; "
        "then "
        "  {params.netMHC}/netMHCIIpan {params.extra} -xlsfile {output} -a {params.alleles} -f {input.peptides} > {log}; "
        "else "
        "  touch {output}; "
        "fi "
        " ) 2> {log}"


rule parse_mhc_out:
    input:
        expand(
            "results/{{mhc}}/{{group}}/{{tumor_alias}}.merged_tumor_normal.{contig}.{{peptide_type}}.xls",
            contig=contigs,
        ),
    output:
        "results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.{peptide_type}.tsv",
    log:
        "logs/parse_mhc_out/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.{peptide_type}.log",
    script:
        "../scripts/group_mhc_output.py"


# rule parse_mhcflurry:
#     input:
#         expand("results/mhcflurry/{{sample}}/{contig}/output.{{peptide_type}}.csv", contig=contigs)
#     output:
#         "results/mhcflurry/{sample}/{sample}.mhc.{peptide_type}.csv"
#     wildcard_constraints:
#         group="wt|mt"
#     log:
#         "logs/parse-mhc/mhcflurry-{sample}-{peptide_type}.log"
#     conda:
#         "../envs/xsv.yaml"
#     shell:
#         "xsv cat rows -d ',' {input} | cut --complement -f2,7,8 > {output}"


rule mhc_csv_table:
    input:
        info="results/microphaser/info/filtered/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.tsv",
        neo="results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.neo.tsv",
        normal="results/{mhc}/{group}.{tumor_alias}.merged_tumor_normal.mhc.normal.tsv",
    output:
        report(
            "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.DNA.tsv",
            caption="../report/WES_results.rst",
            category="Results WES (netMHC)",
        ),
    log:
        "logs/mhc_csv_table/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.log",
    script:
        "../scripts/merge_data.py"


# rule mhcflurry_table:
#     input:
#         info="results/microphaser/info/{sample}/filtered/mhcflurry/{sample}.tsv",
#         neo="results/mhcflurry/{sample}/{sample}.mhc.neo.tsv",
#         normal="results/mhcflurry/{sample}/{sample}.mhc.normal.tsv"
#     output:
#         report("results/neoantigens/mhcflurry/{sample}.WES.tsv", caption="../report/WES_results.rst", category="Results WES (MHCFlurry)")
#     script:
#         "../scripts/merge_mhcflurry.py"


rule add_RNA_info:
    input:
        counts="results/kallisto/{group}.{tumor_alias}",
        table="results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.DNA.tsv",
    output:
        report(
            "results/neoantigens/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.RNA.tsv",
            caption="../report/RNA_results.rst",
            category="Results RNA",
        ),
    params:
        abundance=lambda wc, input: "{}/abundance.tsv".format(input.counts),
    log:
        "logs/add-RNA/{group}.{tumor_alias}.merged_tumor_normal.{mhc}.log",
    script:
        "../scripts/add_rna_info.py"
