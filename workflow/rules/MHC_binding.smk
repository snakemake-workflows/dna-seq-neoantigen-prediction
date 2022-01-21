# rule mhcflurry:
#     input:
#         peptides="results/microphaser/fasta/{sample}/filtered/{sample}.{chr}.{peptide_type}.fa",
#         alleles="results/optitype/{sample}/hla_alleles_{sample}.tsv",
#         wt_alleles=get_germline_optitype
#     output:
#         "results/mhcflurry/{sample}/{chr}/output.{peptide_type}.csv"
#     log:
#         "logs/mhcflurry/{sample}-{chr}-{peptide_type}.log"
#     run:
#         if "wt" in input.peptides:
#             alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
#         else:
#             alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
#         cmd = "if [ -s {input.peptides} ]; then mhctools --mhc-predictor mhcflurry --mhc-alleles {alleles} --input-fasta-file {input.peptides} --output-csv {output} > {log}; else touch {output}; fi"
#         shell(cmd)


rule netMHCpan:
    input:
        peptides="results/microphaser/fasta/{cancer_sample}/filtered/netMHCpan/{cancer_sample}.{chr}.{peptide_type}.fa",
        alleles=get_alleles_MHCI,
    output:
        "results/netMHCpan/{cancer_sample}/{chr}/{cancer_sample}.{chr}.{peptide_type}.xls",
    log:
        "logs/netMHCpan/{cancer_sample}-{chr}-{peptide_type}.log",
    params:
        extra=config["affinity"]["netMHCpan"]["params"],
        netMHC=config["affinity"]["netMHCpan"]["location"],
    run:
        alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
        cmd = "if [ -s {input.peptides} ]; then {params.netMHC}/netMHCpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)


rule netMHCIIpan:
    input:
        peptides="results/microphaser/fasta/{cancer_sample}/filtered/netMHCIIpan/{cancer_sample}.{chr}.{peptide_type}.fa",
        alleles=get_alleles_MHCII,
    output:
        "results/netMHCIIpan/{cancer_sample}/{chr}/{cancer_sample}.{chr}.{peptide_type}.xls",
    log:
        "logs/netMHCIIpan/{cancer_sample}-{chr}-{peptide_type}.log",
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
            "results/{{mhc}}/{{cancer_sample}}/{chr}/{{cancer_sample}}.{chr}.{{peptide_type}}.xls",
            chr=contigs,
        ),
    output:
        "results/{mhc}/{cancer_sample}/{cancer_sample}.mhc.{peptide_type}.tsv",
    log:
        "logs/parse-mhc/{mhc}-{cancer_sample}-{peptide_type}.log",
    wildcard_constraints:
        group="wt|mt",
    script:
        "../scripts/group_mhc_output.py"


# rule parse_mhcflurry:
#     input:
#         expand("results/mhcflurry/{{sample}}/{chr}/output.{{peptide_type}}.csv", chr=contigs)
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
        info="results/microphaser/info/{cancer_sample}/filtered/{mhc}/{cancer_sample}.tsv",
        neo="results/{mhc}/{cancer_sample}/{cancer_sample}.mhc.neo.tsv",
        normal="results/{mhc}/{cancer_sample}/{cancer_sample}.mhc.normal.tsv",
    output:
        report(
            "results/neoantigens/{mhc}/{cancer_sample}.DNA.tsv",
            caption="../report/WES_results.rst",
            category="Results WES (netMHC)",
        ),
    log:
        "logs/create-mhc-table/{mhc}-{cancer_sample}.log",
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
        counts="results/kallisto/{cancer_sample}",
        table="results/neoantigens/{mhc}/{cancer_sample}.DNA.tsv",
    output:
        report(
            "results/neoantigens/{mhc}/{cancer_sample}.RNA.tsv",
            caption="../report/RNA_results.rst",
            category="Results RNA",
        ),
    params:
        abundance=lambda wc, input: "{}/abundance.tsv".format(input.counts),
    log:
        "logs/add-RNA/{mhc}-{cancer_sample}.log",
    script:
        "../scripts/add_rna_info.py"
