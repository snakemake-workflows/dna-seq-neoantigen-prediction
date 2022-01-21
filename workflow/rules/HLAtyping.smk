rule HLA_LA:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai",
        index="resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
    output:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt",
    threads: 7
    log:
        "logs/HLA-LA/{sample}.log",
    params:
        graph=lambda w, input: os.path.basename(os.path.dirname(input.index)),
        graphdir=lambda w, input: os.path.dirname(os.path.dirname(input.index)),
    conda:
        "../envs/hla_la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.sample} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir results/HLA-LA/output --maxThreads {threads} > {log} 2>&1"


rule parse_HLA_LA:
    input:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt",
    output:
        report(
            "results/HLA-LA/hlaI_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
        report(
            "results/HLA-LA/hlaII_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
    log:
        "logs/parse-HLA-LA/{sample}.log",
    script:
        "../scripts/parse_HLA_types.py"


rule razers3:
    input:
        reads="results/merged/DNA/{sample}_{read}.fastq.gz",
    output:
        bam="results/razers3/bam/{sample}_{read}.bam",
    threads: 8
    log:
        "logs/razers3/{sample}_{read}.log",
    params:
        genome=config["HLAtyping"]["optitype_data"],
        extra=config["params"]["razers3"],
    wrapper:
        "0.61.0/bio/razers3"


rule bam2fq:
    input:
        "results/razers3/bam/{sample}_{read}.bam",
    output:
        "results/razers3/fastq/{sample}_{read}.fished.fastq",
    params:
        "",
    log:
        "logs/razers3-bam2fq/{sample}-{read}.log",
    threads: 1
    wrapper:
        "0.61.0/bio/samtools/bam2fq/interleaved"


rule OptiType:
    input:
        reads=get_optitype_reads_input,
    output:
        multiext(
            "results/optitype/{sample}/{sample}", "_coverage_plot.pdf", "_result.tsv"
        ),
    log:
        "logs/optitype/{sample}.log",
    params:
        extra=config["params"]["optitype"],
        sequencing_type="dna",
    wrapper:
        "0.63.0/bio/optitype"


rule parse_Optitype:
    input:
        "results/optitype/{sample}/{sample}_result.tsv",
    output:
        report(
            "results/optitype/{sample}/hla_alleles_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(Optitype)",
        ),
    log:
        "logs/parse-optitype/{sample}.log",
    shell:
        "cut {input} -f2-7 | awk 'NR == 1 {{print}} NR>1 {{for (i = 1; i<=6; ++i) sub(/^/, \"&HLA-\", $i); print}}' "
        '| sed -e s/[*,:]//g | sed "s/ /\t/g" > {output} 2> {log}'
