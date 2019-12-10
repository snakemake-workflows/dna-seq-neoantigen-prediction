rule HLA_LA:
    input:
        bam="bwa/{sample}.rmdup.bam",
        bai="bwa/{sample}.rmdup.bam.bai"
    output:
        "HLA-LA/output/{sample}/hla/R1_bestguess_G.txt"
    threads: 7
    params:
        graph="PRG_MHC_GRCh38_withIMGT",
        graphdir="/vol/tiny/MaMel-Neoantigens/HLA-LA_graphs",
        extra_refs="../HLA-LA_graphs/additionalReferences/PRG_MHC_GRCh38_withIMGT"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.sample} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir HLA-LA/output --maxThreads {threads}"

rule parse_HLA_LA:
    input:
        "HLA-LA/output/{sample}/hla/R1_bestguess_G.txt"
    output:
        report("HLA-LA/hlaI_{sample}.tsv", caption="../report/HLA-LA_Types.rst", category="HLA-Typing(HLA-LA)"),
        report("HLA-LA/hlaII_{sample}.tsv", caption="../report/HLA-LA_Types.rst", category="HLA-Typing(HLA-LA)")
    script:
        "../scripts/parse_HLA_types.py"

rule razers3:
    input:
        get_reads
    output:
        bam="razers3/bam/{sample}_{group}.bam",
        fastq="razers3/fastq/{sample}_{group}.fished.fastq"
    threads: 10
    params:
        hlaref=config["reference"]["hla_data"],
        extra=config["params"]["razers3"]
    conda:
        "../envs/optitype.yaml"
    shell:
        "razers3 -tc {threads} {params.extra} {params.hlaref} {input} -o {output.bam} && samtools bam2fq {output.bam} > {output.fastq}"

rule OptiType:
    input:
        f1='razers3/fastq/{sample}_1.fished.fastq',
        f2='razers3/fastq/{sample}_2.fished.fastq'
    output:
        report("optitype/{sample}/hla_alleles_{sample}.tsv", caption="../report/HLA_Types.rst", category="HLA-Typing")
    params:
        outdir="optitype/{sample}/",
        conf=config["params"]["optitype"]
    conda:
        "../envs/optitype.yaml"
    shell:
        "OptiTypePipeline.py -i {input.f1} {input.f2} --dna --outdir {params.outdir} -c {params.conf} "
        "&& cat {params.outdir}*/*_result.tsv | cut - -f2-7 | awk 'NR == 1 {{print}} NR>1 {{for (i = 1; i<=6; ++i) sub(/^/, \"&HLA-\", $i); print}}' "
        "| sed -e s/[*,:]/''/g | sed s/' '/'\t'/g > {output}"

rule mhcflurry:
    input:
        peptides="microphaser/fasta/{sample}/filtered/{sample}.{chr}.{group}.fa",
        alleles="optitype/{sample}/hla_alleles_{sample}.tsv",
        wt_alleles=get_germline_optitype
    output:
        "mhcflurry/{sample}/{chr}/output.{group}.csv"
    log:
        "logs/mhcflurry/{sample}/{chr}/log.{group}.txt"
#    conda:
#        "../envs/mhctools.yaml"
    run:
        if "wt" in input.peptides:
            alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
        else:
            alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
        cmd = "if [ -s {input.peptides} ]; then mhctools --mhc-predictor mhcflurry --mhc-alleles {alleles} --input-fasta-file {input.peptides} --output-csv {output} > {log}; else touch {output}; fi"
        shell(cmd)

rule netMHCpan:
    input:
        peptides="microphaser/fasta/{sample}/filtered/{sample}.{chr}.{group}.fa",
        alleles="optitype/{sample}/hla_alleles_{sample}.tsv",
        wt_alleles=get_germline_optitype
    output:
        "netMHCpan/{sample}/{chr}/{sample}.{chr}.{group}.xls",
    log:
        "logs/netMHCpan/{sample}/{chr}/{sample}.{chr}.{group}.log"
    params:
        extra = config["params"]["netMHCpan"]
    run:
        if "wt" in input.peptides:
            alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
        else:
            alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
        cmd = "if [ -s {input.peptides} ]; then ../netMHCpan-4.0/netMHCpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)


rule netMHC2:
    input:
        peptides="microphaser/fasta/{sample}/filtered/{sample}.{chr}.{group}.fa",
        alleles = "HLA-LA/hlaII_{sample}.tsv",
        wt_alleles=get_germline_hla
    output:
        "netMHC2pan/{sample}/{chr}/{sample}.{chr}.{group}.xls",
    log:
        "logs/netMHC2pan/{sample}/{chr}/{sample}.{chr}.{group}.log"
    params:
        extra=config["params"]["netMHCIIpan"]
    run:
        if "wt" in input.peptides:
            alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t")["Allele"].tolist())
        else:
            alleles = ",".join(pd.read_csv(input.alleles, sep="\t")["Allele"].tolist())
        cmd = "if [ -s {input.peptides} ]; then ../netMHCIIpan-3.2/netMHCIIpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)


rule parse_mhc_out:
    input:
        expand("{{mhc}}/{{sample}}/{chr}/{{sample}}.{chr}.{{group}}.xls", chr=CHROMOSOMES)
    output:
        "{mhc}/{sample}/{sample}.mhc.{group}.tsv"
    wildcard_constraints:
        group="wt|mt"
    script:
        "../scripts/group_mhc_output.py"

rule parse_mhcflurry:
    input:
        expand("mhcflurry/{{sample}}/{chr}/output.{{group}}.csv", chr=CHROMOSOMES)
    output:
        "mhcflurry/{sample}/{sample}.mhc.{group}.csv"
    wildcard_constraints:
        group="wt|mt"
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d ',' {input} | cut --complement -f2,7,8 > {output}"

rule mhc_csv_table:
    input:
        info="microphaser/info/{sample}/filtered/{sample}.tsv",
        mt="{mhc}/{sample}/{sample}.mhc.mt.tsv",
        wt="{mhc}/{sample}/{sample}.mhc.wt.tsv"
    output:
        report("results/{mhc}/{sample}.WES.tsv", caption="../report/WES_results.rst", category="Results WES")
    script:
        "../scripts/merge_data.py"

rule mhcflurry_table:
    input:
        info="microphaser/info/{sample}/filtered/{sample}.tsv",
        mt="mhcflurry/{sample}/{sample}.mhc.mt.tsv",
        wt="mhcflurry/{sample}/{sample}.mhc.wt.tsv"
    output:
        report("results/mhcflurry/{sample}.WES.tsv", caption="../report/WES_results.rst", category="Results WES")
    script:
        "../scripts/merge_mhcflurry.py"

rule add_RNA_info:
    input:
        counts="transcriptome/kallisto/{sample}/abundance.tsv",
        table="results/{mhc}/{sample}.WES.tsv"
    output:
        report("results/{mhc}/{sample}.RNA.tsv", caption="../report/RNA_results.rst", category="Results RNA")
    script:
        "../scripts/add_rna_info.py"
