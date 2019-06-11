def get_reads(wildcards):
    return get_seperate(**wildcards)#expand("../raw/{sample}_{group}.fastq.gz",

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
        "optitype/{sample}/hla_alleles.tsv"
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
        peptides="microphaser/fasta/{tumor}-{normal}/filtered/{tumor}-{normal}.{chr}.{group}.fa",
        alleles="optitype/{tumor}/hla_alleles.tsv",
        wt_alleles="optitype/{normal}/hla_alleles.tsv"
    output:
        "mhcflurry/{tumor}-{normal}/{chr}/output.{group}.csv"
    log:
        "logs/mhcflurry/{tumor}-{normal}/{chr}/log.{group}.txt"
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
        peptides="microphaser/fasta/{tumor}-{normal}/filtered/{tumor}-{normal}.{chr}.{group}.fa",
        alleles="optitype/{tumor}/hla_alleles.tsv",
        wt_alleles="optitype/{normal}/hla_alleles.tsv"
    output:
        "netMHCpan/{tumor}-{normal}/{chr}/{tumor}-{normal}.{chr}.{group}.xls",
    log:
        "logs/netMHCpan/{tumor}-{normal}/{chr}/{tumor}-{normal}.{chr}.{group}.log"
    params:
        extra = config["params"]["netMHCpan"]
    run:
        if "wt" in input.peptides:
            alleles = ",".join(pd.read_csv(input.wt_alleles, sep="\t").iloc[0])
        else:
            alleles = ",".join(pd.read_csv(input.alleles, sep="\t").iloc[0])
        cmd = "if [ -s {input.peptides} ]; then ../netMHCpan-4.0/netMHCpan {params.extra} -xlsfile {output} -a {alleles} -f {input.peptides} > {log}; else touch {output}; fi"
        shell(cmd)

#rule netMHC2:
#    input:
#        peptides = "out/peptides/relevant/{tumor}-{normal}/{tumor}-{normal}.{C}.filtered.pep",
#        alleles = "out/alleles_MHC2/{tumor}_result.tsv"
#    output:
#        xls = "out/mhc2_binding/relevant/{tumor}-{normal}/{C}/{tumor}-{normal}_{C}.xls",
#        log = "out/mhc2_binding/relevant/{tumor}-{normal}/{C}/{tumor}-{normal}_{C}.log"
#    run:
#        #shell('mkdir out/mhc_binding/{wildcards.P}/{wildcards.C}')
#        alleles=open(input.alleles)
#        head=next(alleles)
#        line=next(alleles)
#        if line.startswith('H'):
#            line = ','.join(line.rstrip().split('\t'))
#            #line = ','.join(list(set(line)))
#            cmd="netMHC2pan -BA -l 9 -s -xls -xlsfile {output.xls} -a " + line + " -f {input.peptides} > {output.log}"
#            shell(cmd)


rule parse_wt_mhc_out:
    input:
        expand("netMHCpan/{{tumor}}-{{normal}}/{chr}/{{tumor}}-{{normal}}.{chr}.{{group}}.xls", chr=CHROMOSOMES)
    output:
        "netMHCpan/{tumor}-{normal}/{tumor}-{normal}.mhc.{group}.tsv"
    wildcard_constraints:
        group="wt|mt"
    run:
        s = input[0]
        shell("python scripts/group_mhc_output.py " + s + " >> {output}")
        for i in input[1:]:
            shell("python scripts/group_mhc_output.py " + i + " | grep -v -h '^Pos' - >> {output}") # xsv

rule parse_mhcflurry:
    input:
        expand("mhcflurry/{{tumor}}-{{normal}}/{chr}/output.{{group}}.csv", chr=CHROMOSOMES)
    output:
        "mhcflurry/{tumor}-{normal}/{tumor}-{normal}.mhc.{group}.csv"
    wildcard_constraints:
        group="wt|mt"
    conda:
        "../envs/xsv.yaml"
    shell:
        "xsv cat rows -d ',' {input} | cut --complement -f2,7,8 > {output}"

rule mhc_csv_table:
    input:
        info="microphaser/info/{tumor}-{normal}/filtered/{tumor}-{normal}.tsv",
        mt="netMHCpan/{tumor}-{normal}/{tumor}-{normal}.mhc.mt.tsv",
        wt="netMHCpan/{tumor}-{normal}/{tumor}-{normal}.mhc.wt.tsv"
    output:
        "results/netMHCpan/{tumor}-{normal}.tsv"
    script:
        "../scripts/merge_data.py"

rule mhcflurry_table:
    input:
        info="microphaser/info/{tumor}-{normal}/filtered/{tumor}-{normal}.tsv",
        mt="mhcflurry/{tumor}-{normal}/{tumor}-{normal}.mhc.mt.csv",
        wt="mhcflurry/{tumor}-{normal}/{tumor}-{normal}.mhc.wt.csv"
    output:
        "results/mhcflurry/{tumor}-{normal}.tsv"
    script:
        "../scripts/merge_mhcflurry.py"
