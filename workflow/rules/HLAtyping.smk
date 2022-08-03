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


rule net_mhc_pan_alleles:
    output:
        mhc_one_alleles="resources/hla_alleles/available_alleles.net_mhc_pan.txt",
    conda:
        "../envs/tcsh.yaml"
    log:
        "logs/net_mhc_pan/available_alleles.net_mhc_pan.log",
    params:
        net_mhc=config["params"]["net_mhc_pan"]["location"],
    shell:
        "{params.net_mhc}/net_mhc_pan -listMHC > {output.mhc_one_alleles} 2> {log}"


rule net_mhc_two_pan_alleles:
    output:
        mhc_two_alleles="resources/hla_alleles/available_alleles.net_mhc_two_pan.txt",
    conda:
        "../envs/tcsh.yaml"
    log:
        "logs/net_mhc_pan/available_alleles.net_mhc_two_pan.log",
    params:
        net_mhc=config["params"]["net_mhc_two_pan"]["location"],
    shell:
        "{params.net_mhc}/net_mhc_two_pan -list > {output.mhc_two_alleles} 2> {log}"


rule parse_and_filter_hla_alleles_for_netmhc:
    input:
        hla_la_bestguess="results/HLA-LA/output/{group}_{alias}/hla/R1_bestguess_G.txt",
        mhc_one_alleles="resources/hla_alleles/available_alleles.net_mhc_pan.txt",
        mhc_two_alleles="resources/hla_alleles/available_alleles.net_mhc_two_pan.txt",
    output:
        hlaI=report(
            "results/HLA-LA/{group}.{alias}.hlaI.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
        hlaII=report(
            "results/HLA-LA/{group}.{alias}.hlaII.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
    log:
        "logs/parse-HLA-LA/{group}.{alias}.log",
    script:
        "../scripts/parse_and_filter_hla_alleles_for_netmhc.py"
