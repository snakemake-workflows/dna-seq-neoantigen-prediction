rule strelka_tumor:
    input:
        normal=get_normal_bam(),
        normal_index=get_normal_bam(ext=".bam.bai"),
        tumor=get_tumor_bam(),
        tumor_index=get_tumor_bam(ext=".bam.bai"),
        fasta="resources/genome.fasta",
        fasta_index="resources/genome.fasta.fai",
        callregions="resources/genome.callregions.bed.gz",
    output:
        "results/strelka/{group}.strelka_somatic.snvs.vcf.gz",
        "results/strelka/{group}.strelka_somatic.indels.vcf.gz",
    log:
        "logs/calling/strelka/{group}.strelka_somatic.log",
    params:
        config_extra="--callRegions {} {}".format(
            "resources/genome.callregions.bed.gz",
            config["params"]["strelka"]["config"],
        ),
        run_extra=config["params"]["strelka"]["run"],
    threads: 22
    wrapper:
        "0.65.0/bio/strelka/somatic"


rule strelka_germline:
    input:
        bam=get_normal_bam(),
        normal_index=get_normal_bam(ext=".bam.bai"),
        fasta="resources/genome.fasta",
        fasta_index="resources/genome.fasta.fai",
        callregions="resources/genome.callregions.bed.gz",
    output:
        "results/strelka/{group}.strelka_germline.variants.vcf.gz",
    log:
        "logs/calling/strelka_germline/{group}.log",
    params:
        config_extra="--callRegions {} {}".format(
            "resources/genome.callregions.bed.gz",
            config["params"]["strelka"]["config"],
        ),
        run_extra="",
    threads: 22
    wrapper:
        "0.65.0/bio/strelka/germline"


rule vcf_to_bcf:
    input:
        "{variants}.vcf.gz",
    output:
        "{variants}.output.bcf",
    log:
        "logs/bcftools/to-bcf/{variants}.log",
    params:
        "-O b -f PASS",
    wrapper:
        "0.60.0/bio/bcftools/view"


rule concat_somatic:
    input:
        calls=expand(
            "results/strelka/{{group}}.strelka_somatic.{type}.output.bcf",
            type=["snvs", "indels"],
        ),
        indices=expand(
            "results/strelka/{{group}}.strelka_somatic.{type}.output.bcf.csi",
            type=["snvs", "indels"],
        ),
    output:
        "results/strelka/{group}.strelka_somatic.bcf",
    log:
        "bcftools/concat_somatic/{group}.log",
    params:
        "-O b -a",
    wrapper:
        "0.60.0/bio/bcftools/concat"


rule get_tumor_from_somatic:
    input:
        "results/strelka/{group}.strelka_somatic.bcf",
    output:
        "results/strelka/{group}.strelka_somatic.tumor.bcf",
    log:
        "logs/bcftools/get_tumor_from_somatic/{group}.strelka_somatic.tumor.log",
    params:
        "-O b -s TUMOR",
    wrapper:
        "0.60.0/bio/bcftools/view"


rule reheader_germline:
    input:
        vcf="results/strelka/{group}.strelka_germline.variants.output.bcf",
        samples="resources/sampleheader.txt",
    output:
        "results/strelka/{group}.strelka_germline.variants.reheader.bcf",
    log:
        "logs/bcftools/reheader_germline/{group}.log",
    params:
        extra="",
        view_extra="-O b",
    wrapper:
        "0.60.0/bio/bcftools/reheader"


rule concat_variants:
    input:
        calls=[ 
            "results/strelka/{group}.strelka_somatic.tumor.bcf",
            "results/strelka/{group}.strelka_germline.variants.reheader.bcf",
        ],
        index=[ 
            "results/strelka/{group}.strelka_somatic.tumor.bcf.csi",
            "results/strelka/{group}.strelka_germline.variants.reheader.bcf.csi",
        ],
    output:
        "results/strelka/merged/{group}.strelka_somatic.strelka_germline.bcf",
    log:
        "bcftools/concat_variants/{group}.strelka_somatic.strelka_germline.log",
    params:
        extra="-O b -a",
    wrapper:
        "0.64.0/bio/bcftools/concat"


rule norm_vcf:
    input:
        "{prefix}.bcf",
        genome="resources/genome.fasta",
    output:
        "{prefix}.norm.bcf",
    log:
        "logs/bcftools/norm/{prefix}.log",
    params:
        "-f {} -O b -m-".format("resources/genome.fasta"),  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.65.0/bio/bcftools/norm"
