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


rule kallisto_index:
    input:
        "resources/genome.cdna.fasta",
    output:
        "resources/kallisto/transcripts.idx",
    params:
        extra="",
    log:
        "logs/kallisto/index.log",
    cache: True
    wrapper:
        "0.60.1/bio/kallisto/index"


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


rule STAR_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        directory("resources/STAR_index"),
    params:
        sjdb_overhang="100",
        extra="",
    log:
        "logs/star/index.log",
    threads: 32
    cache: True
    wrapper:
        "0.42.0/bio/star/index"


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


rule get_callregions:
    input:
        "resources/genome.fasta.fai",
    output:
        "resources/genome.callregions.bed.gz",
    log:
        "logs/get-callregions.log",
    params:
        n_contigs=config["ref"]["n_chromosomes"],
    conda:
        "../envs/htslib.yaml"
    shell:
        "paste <(cut -f1 {input}) <(yes 0 | head -n {params.n_contigs}) <(cut -f2 {input})"
        " | head -n {params.n_contigs} | bgzip -c > {output} && tabix -p bed {output}"


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    cache: True
    wrapper:
        "0.45.1/bio/bwa/index"


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


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    cache: True
    wrapper:
        "0.59.2/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=config["ref"]["release"],
    log:
        "logs/vep/plugins.log",
    cache: True
    wrapper:
        "0.59.2/bio/vep/plugins"


rule make_sampleheader:
    output:
        "resources/sampleheader.txt",
    log:
        "logs/germline-reheader-sample.log",
    shell:
        "echo 'TUMOR' > {output}"
