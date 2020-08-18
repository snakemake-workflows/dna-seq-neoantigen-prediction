rule strelka_somatic:
    input:
        normal=get_normal_bam,
        normal_index=get_normal_bai,
        tumor="results/recal/{sample}.sorted.bam",
        tumor_index="results/recal/{sample}.sorted.bam.bai",
        fasta="resources/genome.fasta",
        fasta_index="resources/genome.fasta.fai",
        callregions="resources/genome.callregions.bed.gz"
    output:
        "results/strelka/somatic/{sample}/results/variants/somatic.snvs.vcf.gz",
        "results/strelka/somatic/{sample}/results/variants/somatic.indels.vcf.gz"
    log:
        "logs/calling/strelka_somatic/{sample}.log"
    conda:
        "../envs/strelka.yaml"
    params:
        config_extra="--callRegions resources/genome.callregions.bed.gz {}".format(config["params"]["strelka"]),
        run_extra=""
    threads: 22
    shell:
        "(configureStrelkaSomaticWorkflow.py "  # Configuration script
        "--normalBam {input.normal} "  # Path to normal bam (if any)
        "--tumorBam {input.tumor} "  # Path to tumor bam
        "--referenceFasta {input.fasta} "  # Path to fasta file
        "--runDir results/strelka/somatic/{wildcards.sample} "  # Path to output directory
        "{params.config_extra} "  # Extra parametersfor configuration
        " && "
        "results/strelka/somatic/{wildcards.sample}/runWorkflow.py "  # Run the pipeline
        "--mode local "  # Stop internal job submission
        "--jobs {threads} "  # Nomber of threads
        "{params.run_extra}) "  # Extra parameters for runWorkflow
        "> {log} 2>&1"  # Logging behaviour"
    # wrapper:
        # "f834d6c44eaa4b534a0467cf1620f998ad289667/bio/strelka/somatic"

rule strelka_germline:
    input:
        bam="results/recal/{normal}.sorted.bam",
        normal_index="results/recal/{normal}.sorted.bam.bai",
        fasta="resources/genome.fasta",
        fasta_index="resources/genome.fasta.fai",
        callregions="resources/genome.callregions.bed.gz"
    output:
        "results/strelka/germline/{normal}/results/variants/variants.vcf.gz"
    log:
        "logs/calling/strelka_germline/{normal}.log"
    params:
        config_extra="--callRegions resources/genome.callregions.bed.gz {}".format(config["params"]["strelka"]),
        run_extra=""
    threads: 22
    wrapper:
        "0.60.0/bio/strelka/germline"

rule vcf_to_bcf:
    input:
        "{variants}.vcf.gz"
    output:
        "{variants}.output.bcf"
    params:
        "-O b -f PASS"
    wrapper:
        "0.60.0/bio/bcftools/view"

rule index_bcf:
    input:
        "{variants}.bcf"
    output:
        "{variants}.bcf.csi"
    params:
        extra=""
    wrapper:
        "0.60.0/bio/bcftools/index"


rule concat_somatic:
    input:
        calls=expand("results/strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf", type=["snvs","indels"]),
        indices=expand("results/strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf.csi", type=["snvs","indels"])
    output:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.bcf"
    params:
        "-O b -a"
    wrapper:
        "0.60.0/bio/bcftools/concat"

rule get_tumor_from_somatic:
    input:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.bcf"
    output:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf"
    params:
        "-O b -s TUMOR"
    wrapper:
        "0.60.0/bio/bcftools/view"

rule clean_germline:
    input:
        vcf="{germline}/variants.output.bcf",
        samples="config/newsamples.txt"
    output:
        "{germline}/variants.reheader.bcf"
    params:
        extra="",
        view_extra="-O b"
    wrapper:
        "0.60.0/bio/bcftools/reheader"

rule concat_variants:
    input:
        germline=get_germline_variants,
        index_g=get_germline_variants_index,
        somatic="results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf",
        index_s="results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf.csi",
        db="resources/snpEff/GRCh38.86"
    output:
        "results/strelka/merged/{sample}/all_variants.bcf"
    params:
        data_dir=lambda _, input: Path(input.db).parent.resolve(),
        extra="-O v"
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools concat {params.extra} -a {input.somatic} {input.germline} | snpEff -Xmx4g -nodownload -dataDir {params.data_dir} GRCh38.86 - | bcftools view -O u - > {output}"

rule preprocess_variants:
    input:
        variants="{variants}.bcf"
    output:
        "{variants}.prepy.bcf"
    params:
        extra="-L --somatic",
        genome="resources/genome.fasta",
    threads: 2
    wrapper:
        "0.60.0/bio/hap.py/pre.py"
