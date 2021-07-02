# rule mhcflurry:
#     input:
#         peptides="results/microphaser/fasta/{sample}/filtered/{sample}.{chr}.{group}.fa",
#         alleles="results/optitype/{sample}/hla_alleles_{sample}.tsv",
#         wt_alleles=get_germline_optitype
#     output:
#         "results/mhcflurry/{sample}/{chr}/output.{group}.csv"
#     log:
#         "logs/mhcflurry/{sample}-{chr}-{group}.log"
#     run:
#         if "wt" in input.peptides:
#             alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
#         else:
#             alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
#         cmd = "if [ -s {input.peptides} ]; then mhctools --mhc-predictor mhcflurry --mhc-alleles {alleles} --input-fasta-file {input.peptides} --output-csv {output} > {log}; else touch {output}; fi"
#         shell(cmd)


rule netMHCpan:
    input:
        peptides="results/microphaser/fasta/{sample}/filtered/netMHCpan/{sample}.{chr}.{group}.fa",
        alleles=get_alleles_MHCI,
    output:
        "results/netMHCpan/{sample}/{chr}/{sample}.{chr}.{group}.xls",
    log:
        "logs/netMHCpan/{sample}-{chr}-{group}.log",
    params:
        extra=config["affinity"]["netMHCpan"]["params"],
        netMHC=config["affinity"]["netMHCpan"]["location"],
    run:
        alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
        cmd = "if [ -s {input.peptides} ]; then {params.netMHC}/netMHCpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)


rule netMHCIIpan:
    input:
        peptides="results/microphaser/fasta/{sample}/filtered/netMHCIIpan/{sample}.{chr}.{group}.fa",
        alleles=get_alleles_MHCII,
    output:
        "results/netMHCIIpan/{sample}/{chr}/{sample}.{chr}.{group}.xls",
    log:
        "logs/netMHCIIpan/{sample}-{chr}-{group}.log",
    params:
        extra=config["affinity"]["netMHCIIpan"]["params"],
        netMHC=config["affinity"]["netMHCIIpan"]["location"],
    run:
        alleles = ",".join(pd.read_csv(input.alleles, sep="\t")["Allele"].tolist())
        cmd = "if [ -s {input.peptides} ]; then {params.netMHC}/netMHCIIpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)


rule parse_mhc_out:
    input:
        expand(
            "results/{{mhc}}/{{sample}}/{chr}/{{sample}}.{chr}.{{group}}.xls",
            chr=contigs,
        ),
    output:
        "results/{mhc}/{sample}/{sample}.mhc.{group}.tsv",
    log:
        "logs/parse-mhc/{mhc}-{sample}-{group}.log",
    wildcard_constraints:
        group="wt|mt",
    script:
        "../scripts/group_mhc_output.py"


# rule parse_mhcflurry:
#     input:
#         expand("results/mhcflurry/{{sample}}/{chr}/output.{{group}}.csv", chr=contigs)
#     output:
#         "results/mhcflurry/{sample}/{sample}.mhc.{group}.csv"
#     wildcard_constraints:
#         group="wt|mt"
#     log:
#         "logs/parse-mhc/mhcflurry-{sample}-{group}.log"
#     conda:
#         "../envs/xsv.yaml"
#     shell:
#         "xsv cat rows -d ',' {input} | cut --complement -f2,7,8 > {output}"


rule mhc_csv_table:
    input:
        info="results/microphaser/info/{sample}/filtered/{mhc}/{sample}.tsv",
        mt="results/{mhc}/{sample}/{sample}.mhc.mt.tsv",
        wt="results/{mhc}/{sample}/{sample}.mhc.wt.tsv",
    output:
        report(
            "results/neoantigens/{mhc}/{sample}.DNA.tsv",
            caption="../report/WES_results.rst",
            category="Results WES (netMHC)",
        ),
    log:
        "logs/create-mhc-table/{mhc}-{sample}.log",
    script:
        "../scripts/merge_data.py"


# rule mhcflurry_table:
#     input:
#         info="results/microphaser/info/{sample}/filtered/mhcflurry/{sample}.tsv",
#         mt="results/mhcflurry/{sample}/{sample}.mhc.mt.tsv",
#         wt="results/mhcflurry/{sample}/{sample}.mhc.wt.tsv"
#     output:
#         report("results/neoantigens/mhcflurry/{sample}.WES.tsv", caption="../report/WES_results.rst", category="Results WES (MHCFlurry)")
#     script:
#         "../scripts/merge_mhcflurry.py"


rule add_RNA_info:
    input:
        counts="results/kallisto/{sample}",
        table="results/neoantigens/{mhc}/{sample}.DNA.tsv",
    output:
        report(
            "results/neoantigens/{mhc}/{sample}.RNA.tsv",
            caption="../report/RNA_results.rst",
            category="Results RNA",
        ),
    params:
        abundance=lambda wc, input: "{}/abundance.tsv".format(input.counts),
    log:
        "logs/add-RNA/{mhc}-{sample}.log",
    script:
        "../scripts/add_rna_info.py"
