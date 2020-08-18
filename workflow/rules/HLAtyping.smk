rule HLA_LA:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai"
    output:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt"
    threads: 7
    params:
        graph="PRG_MHC_GRCh38_withIMGT",
        graphdir=config["reference"]["HLA_LA_graphs"],
        #extra_refs="/vol/tiny/MaMel-Neoantigens/HLA-LA_graphs/additionalReferences/PRG_MHC_GRCh38_withIMGT"
    conda:
        "../envs/hla_la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.sample} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir results/HLA-LA/output --maxThreads {threads}"

rule parse_HLA_LA:
    input:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt"
    output:
        report("results/HLA-LA/hlaI_{sample}.tsv", caption="../report/HLA_Types.rst", category="HLA-Typing(HLA-LA)"),
        report("results/HLA-LA/hlaII_{sample}.tsv", caption="../report/HLA_Types.rst", category="HLA-Typing(HLA-LA)")
    script:
        "../scripts/parse_HLA_types.py"

rule razers3:
    input:
        reads=get_reads
    output:
        bam="results/razers3/bam/{sample}_{group}.bam"
    threads: 8
    params:
        genome=config["reference"]["hla_data"],
        extra=config["params"]["razers3"]
    wrapper:
        "0.61.0/bio/razers3"

rule bam2fq:
    input:
        "results/razers3/bam/{sample}_{group}.bam"
    output:
        "results/razers3/fastq/{sample}_{group}.fished.fastq"
    params:
        ""
    threads: 1
    wrapper:
        "0.61.0/bio/samtools/bam2fq/interleaved"

rule OptiType:
    input:
        reads=expand("results/razers3/fastq/{{sample}}_{fq}.fished.fastq", fq=[1,2])
    output:
        multiext("results/optitype/{sample}/{sample}", "_coverage_plot.pdf", "_result.tsv")
    log:
        "logs/optitype/{sample}.log"
    params:
        extra=config["params"]["optitype"],
        sequencing_type="dna"
    wrapper:
        "0.63.0/bio/optitype"
        

rule parse_Optitype:
    input:
        "results/optitype/{sample}/{sample}_result.tsv"
    output:
        report("results/optitype/{sample}/hla_alleles_{sample}.tsv", caption="../report/HLA_Types.rst", category="HLA-Typing(Optitype)")
    shell:
        "cut {input} -f2-7 | awk 'NR == 1 {{print}} NR>1 {{for (i = 1; i<=6; ++i) sub(/^/, \"&HLA-\", $i); print}}' "
        "| sed -e s/[*,:]/''/g | sed s/' '/'\t'/g > {output}"
        
