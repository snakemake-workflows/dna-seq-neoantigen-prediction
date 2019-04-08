rule razers3:
    input:
        "../raw/{sample}_{group}.fastq.gz"
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
        configfile=config["params"]["optitype"]
    conda:
        "../envs/optitype.yaml"
    shell:
        "OptiTypePipeline.py -i {input.f1} {input.f2} --dna --outdir {params.outdir} -c {params.configfile} "
        "&& cat {params.outdir}*/*_result.tsv | cut - -f2-7 | awk 'NR == 1 {{print}} NR>1 {{for (i = 1; i<=6; ++i) sub(/^/, \"&HLA-\", $i); print}}' "
        "| sed -e s/[*,:]/''/g | sed s/' '/'\t'/g > {output}"



rule netMHCpan:
    input:
        peptides="microphaser/fasta/{tumor}_{normal}/filtered/{tumor}_{normal}.{chr}.{group}.fa",
        alleles="optitype/{tumor}/hla_alleles.tsv",
        wt_alleles="optitype/{normal}/hla_alleles.tsv"
    output:
        "netMHCpan/{tumor}_{normal}/{chr}/{tumor}_{normal}.{chr}.{group}.xls",
    log:
        "logs/netMHCpan/{tumor}_{normal}/{chr}/{tumor}_{normal}.{chr}.{group}.log"
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
#        peptides = "out/peptides/relevant/{tumor}_{normal}/{tumor}_{normal}.{C}.filtered.pep",
#        alleles = "out/alleles_MHC2/{tumor}_result.tsv"
#    output:
#        xls = "out/mhc2_binding/relevant/{tumor}_{normal}/{C}/{tumor}_{normal}_{C}.xls",
#        log = "out/mhc2_binding/relevant/{tumor}_{normal}/{C}/{tumor}_{normal}_{C}.log"
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
        expand("netMHCpan/{{tumor}}_{{normal}}/{chr}/{{tumor}}_{{normal}}.{chr}.{{type}}.xls", chr=CHROMOSOMES)
    output:
        "netMHCpan/{tumor}_{normal}/{tumor}_{normal}.mhc.{type}.tsv"
    wildcard_constraints:
        type="wt|mt"
    run:
        s = input[0]
        shell("python scripts/group_mhc_output.py " + s + " >> {output}")
        for i in input[1:]:
            shell("python scripts/group_mhc_output.py " + i + " | grep -v -h '^Pos' - >> {output}") ## xsv

rule mhc_csv_table:
    input:
        info="microphaser/info/{tumor}_{normal}/filtered/{tumor}_{normal}.tsv",
        mt="netMHCpan/{tumor}_{normal}/{tumor}_{normal}.mhc.mt.tsv",
        wt="netMHCpan/{tumor}_{normal}/{tumor}_{normal}.mhc.wt.tsv"
    output:
        "tables/{tumor}_{normal}.tsv"
    script:
        "../scripts/merge_data.py"
