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


# TODO: remove this rule, once microphaser is fixed to make gene_name optional
rule remove_records_with_gene_name_missing:
    input:
        "resources/genome.gtf",
    output:
        "resources/genome.records_with_gene_name.gtf",
    log:
        "logs/remove_records_with_gene_name_missing.log",
    shell:
        '( grep "gene_name" {input} > {output} ) 2> {log}'


rule split_annotation:
    input:
        "resources/genome.records_with_gene_name.gtf",
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


rule create_somatic_flag_header_line:
    output:
        "resources/somatic_flag_header_line.txt",
    log:
        "logs/create_somatic_flag_header_line.log"
    shell:
        """
        ( echo '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic tumor variant">' > {output} ) 2> {log}
        """


rule create_genome_somatic_flag_bed:
    input:
        "resources/genome.fasta.fai",
    output:
        "resources/genome.somatic_flag.bed",
    log:
        "logs/create_genome_somatic_flag_bed.log"
    conda:
        "../envs/gawk.yaml"
    cache: True
    shell:
        """
        ( awk 'BEGIN {{ OFS="\\t" }} {{ print $1,0,$2 }}' {input} > {output} ) 2> {log}
        """


rule bgzip_genome_somatic_flag_bed:
    input:
        "resources/genome.somatic_flag.bed",
    output:
        "resources/genome.somatic_flag.bed.gz",
    log:
        "logs/bgzip/genome.somatic_flag.log",
    wrapper:
        "v1.7.0/bio/bgzip"


rule tabix_genome_somatic_flag_bed:
    input:
        "resources/genome.somatic_flag.bed.gz",
    output:
        "resources/genome.somatic_flag.bed.gz.tbi",
    conda:
        "../envs/htslib.yaml"
    log:
        "logs/tabix/genome.somatic_flag.log",
    shell:
        "( tabix -p bed {input} ) 2> {log}"


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
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/knownReferences"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping_PRGonly"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/referenceGenomeSimulations"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/sampledReferenceGenomes"),
        directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/translation"),
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
    params:
        graphs_dir=lambda w, output: output[0].replace(
            "/PRG_MHC_GRCh38_withIMGT/PRG", ""
        ),
    log:
        "logs/download-HLA-LA-graph.log",
    shell:
        "( cd {params.graphs_dir} && wget  http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz "
        "&& tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz && rm  PRG_MHC_GRCh38_withIMGT.tar.gz ) 2> {log}"


rule index_HLALA:
    input:
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
    output:
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH_preGapPathIndex",
    conda:
        "../envs/hla_la.yaml"
    params:
        path=lambda wc, input: os.path.dirname(os.path.dirname(input[0])),
        graph=lambda wc, input: os.path.basename(os.path.dirname(input[0])),
    log:
        "logs/index-HLA-LA-graph.log",
    shell:
        "HLA-LA.pl --prepareGraph 1 --customGraphDir {params.path} --graph {params.graph} > {log} 2>&1"


rule index_HLALA_extended_ref:
    input:
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa",
    output:
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.amb",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.ann",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.bwt",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.pac",
        "resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.sa",
    conda:
        "../envs/hla_la.yaml"
    log:
        "logs/index_HLA-LA_extended_ref.log",
    shell:
        "bwa index {input} > {log} 2>&1"


rule make_sampleheader:
    output:
        "resources/sampleheader.txt",
    log:
        "logs/germline-reheader-sample.log",
    shell:
        "echo 'TUMOR' > {output}"
