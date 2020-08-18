rule get_genome:
    output:
        "resources/genome.fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_cdna:
    output:
        "resources/genome.cdna.fasta"
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf"
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    cache: True
    log:
        "logs/get_annotation.log"
    wrapper:
        "0.45.1/bio/reference/ensembl-annotation"


rule split_annotation:
    input:
        "resources/genome.gtf"
    output:
        "resources/annotation/{contig}.gtf"
    shell:
       "grep '^{wildcards.contig}' {input} > {output}"
        #"awk '!/^#/{{print >\"{output}/\"$1\".gtf\"}}' {input}"


rule genome_faidx:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    cache: True
    wrapper:
        "0.45.1/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.dict"
    log:
        "logs/picard/create_dict.log"
    cache: True
    wrapper:
        "0.45.1/bio/picard/createsequencedictionary"

rule get_callregions:
    input:
        "resources/genome.fasta.fai",
    output:
        "resources/genome.callregions.bed.gz"
    params:
        n_contigs=config["ref"]["n_chromosomes"]
    conda:
        "../envs/index.yaml"
    shell:
        "paste <(cut -f1 {input}) <(yes 0 | head -n {params.n_contigs}) <(cut -f2 {input})"
        " | head -n {params.n_contigs} | bgzip -c > {output} && tabix -p bed {output}"


#rule get_known_variants:
#    input:
#        # use fai to annotate contig lengths for GATK BQSR
#        fai="resources/genome.fasta.fai"
#    output:
#        vcf="resources/variation.vcf.gz"
#    params:
#        species=config["ref"]["species"],
#        release=config["ref"]["release"],
#        type="all"
#    cache: True
#    wrapper:
#        "0.45.1/bio/reference/ensembl-variation"


#rule remove_iupac_codes:
#    input:
#        "resources/variation.vcf.gz"
#    output:
#        "resources/variation.noiupac.vcf.gz"
#    conda:
#        "../envs/rbt.yaml"
#    cache: True
#    shell:
#        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


#rule tabix_known_variants:
#    input:
#        "resources/{prefix}.vcf.gz"
#    output:
#        "resources/{prefix}.vcf.gz.tbi"
#    params:
#        "-p vcf"
#    log:
#        "logs/tabix/{prefix}.log"
#    cache: True
#    wrapper:
#        "0.45.1/bio/tabix"


rule bwa_index:
    input:
        "resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    cache: True
    wrapper:
        "0.45.1/bio/bwa/index"


rule get_snpeff_data:
    output:
        directory("resources/snpEff/{reference}")
    log:
        "logs/snpeff/download/{reference}.log"
    cache: True
    params:
        reference="{reference}"
    wrapper:
        "0.60.1/bio/snpeff/download"

rule get_vep_cache:
    output:
        directory("resources/vep/cache")
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    log:
        "logs/vep/cache.log"
    wrapper:
        "0.59.2/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=config["ref"]["release"]
    log:
        "logs/vep/plugins.log"
    wrapper:
        "0.59.2/bio/vep/plugins"
