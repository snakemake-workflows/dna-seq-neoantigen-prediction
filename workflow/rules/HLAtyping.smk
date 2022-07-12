rule HLA_LA:
    input:
        bam=get_bam_from_group_and_alias(),
        bai=get_bam_from_group_and_alias(ext=".bai"),
        index="resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
        ext_idx="resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.pac",
    output:
        "results/HLA-LA/output/{group}_{alias}/hla/R1_bestguess_G.txt",
    threads: 7
    log:
        "logs/HLA-LA/{group}_{alias}.log",
    params:
        graph=lambda w, input: os.path.basename(os.path.dirname(input.index)),
        graphdir=lambda w, input: os.path.dirname(os.path.dirname(input.index)),
        workdir=lambda w, output: os.path.dirname(
            os.path.dirname(os.path.dirname(output[0]))
        ),
    conda:
        "../envs/hla_la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.group}_{wildcards.alias} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir {params.workdir} --maxThreads {threads} > {log} 2>&1"


rule parse_HLA_LA:
    input:
        "results/HLA-LA/output/{group}_{alias}/hla/R1_bestguess_G.txt",
    output:
        report(
            "results/HLA-LA/{group}.{alias}.hlaI.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
        report(
            "results/HLA-LA/{group}.{alias}.hlaII.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
    log:
        "logs/parse-HLA-LA/{group}.{alias}.log",
    script:
        "../scripts/parse_HLA_types.py"
