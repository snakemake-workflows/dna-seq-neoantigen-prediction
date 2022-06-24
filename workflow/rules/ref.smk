rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_cdna:
    output:
        "resources/genome.cdna.fasta",
    log:
        "logs/get-cdna.log",
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    cache: True
    log:
        "logs/get-annotation.log",
    wrapper:
        "0.45.1/bio/reference/ensembl-annotation"


rule split_annotation:
    input:
        "resources/genome.gtf",
    output:
        "resources/annotation/{contig}.gtf",
    log:
        "logs/split-annotation.{contig}.log",
    shell:
        "grep '^{wildcards.contig}\t' {input} > {output}"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.45.1/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/picard/create-dict.log",
    cache: True
    wrapper:
        "0.45.1/bio/picard/createsequencedictionary"


rule download_HLALA_graph:
    output:
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/PRG"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/knownReferences"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping_PRGonly"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/referenceGenomeSimulations"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/sampledReferenceGenomes"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/translation"),
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
    log:
        "logs/download-HLA-LA-graph.log",
    shell:
        "cd resources/graphs && wget  http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz "
        "&& tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz"


rule index_HLALA:
    input:
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
    output:
        multiext(
            "resources/graphs/PRG_MHC_GRCh38_withIMGT/",
            "serializedGRAPH",
            "serializedGRAPH_preGapPathindex",
        ),
    cache: True
    conda:
        "../envs/hla_la.yaml"
    params:
        path=lambda wc, input: os.path.dirname(os.path.dirname(input[0])),
        graph=lambda wc, input: os.path.basename(os.path.dirname(input[0])),
    log:
        "logs/index-HLA-LA-graph.log",
    shell:
        "HLA-LA.pl --prepareGraph 1 --customGraphDir {params.path} --graph {params.graph} > {log} 2>&1"


rule download_hla_allele_list:
    input:
        HTTP.remote(
            expand(
                "raw.githubusercontent.com/ANHIG/IMGTHLA/{version}/Allelelist.txt",
                version=config["HLAtyping"]["imgt_hla_version"]
            ),
        ),
    output:
        "resources/hla_alleles/Allelelist.txt",
    log:
        "logs/hla_alleles/download_Allelelist.log",
    shell:
        "( mv {input} {output} ) 2> {log}"


rule get_hla_allele_names:
    input:
        "resources/hla_alleles/Allelelist.txt",
    output:
        "resources/hla_alleles/hla_allele_names.txt",
    log:
        "logs/hla_alleles/hla_allele_names.log",
    conda:
        "../envs/grep_sed.yaml"
    shell:
        '( grep -v "^\\(#\\|Allele\\)" {input} | '
        '  cut -d "," -f 2,2 | '
        '  cut -d "*" -f 1,1 | '
        "  uniq | "
        "  sed -e 's/^\\([A-Z]\\)$/HLA-\\1/' | "
        "  sed -e 's/^\\(D[A-Z]\\{{2,2\\}}[1-9]*\\)$/HLA-\\1/' "
        "  >{output} ) 2> {log}"


rule get_hla_regions_from_gtf:
    input:
        gtf="resources/genome.gtf",
        allele_names="resources/hla_alleles/hla_allele_names.txt",
    output:
        "resources/hla_alleles/hla_allele_regions.bed",
    log:
        "logs/hla_alleles/hla_allele_regions.log",
    conda:
        "../envs/rust.yaml"
    script:
        "../scripts/hla_regions_from_gtf.rs"


rule expand_hla_regions:
    input:
        bed="resources/hla_alleles/hla_allele_regions.bed",
        genome="resources/genome.fasta.fai",
    output:
        "resources/hla_alleles/hla_allele_regions.expanded_1000.bed",
    log:
        "logs/hla_alleles/hla_allele_regions.expanded_1000.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "( sort {input.bed} | bedtools slop -b 1000 -g {input.genome} | bedtools merge > {output} ) 2> {log}"


rule make_sampleheader:
    output:
        "resources/sampleheader.txt",
    log:
        "logs/germline-reheader-sample.log",
    shell:
        "echo 'TUMOR' > {output}"
