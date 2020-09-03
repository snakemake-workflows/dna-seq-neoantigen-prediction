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
    # conda:
    #     "../envs/strelka.yaml"
    params:
        config_extra="--callRegions {} {}".format("resources/genome.callregions.bed.gz", config["params"]["strelka"]["config"]),
        run_extra=config["params"]["strelka"]["run"]
    threads: 22
    #shell:
        # "(configureStrelkaSomaticWorkflow.py "  # Configuration script
        # "--normalBam {input.normal} "  # Path to normal bam (if any)
        # "--tumorBam {input.tumor} "  # Path to tumor bam
        # "--referenceFasta {input.fasta} "  # Path to fasta file
        # "--runDir results/strelka/somatic/{wildcards.sample} "  # Path to output directory
        # "{params.config_extra} "  # Extra parametersfor configuration
        # " && "
        # "results/strelka/somatic/{wildcards.sample}/runWorkflow.py "  # Run the pipeline
        # "--mode local "  # Stop internal job submission
        # "--jobs {threads} "  # Nomber of threads
        # "{params.run_extra}) "  # Extra parameters for runWorkflow
        # "> {log} 2>&1"  # Logging behaviour"
    wrapper:
        "0.65.0/bio/strelka/somatic"

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
        config_extra="--callRegions {} {}".format("resources/genome.callregions.bed.gz", config["params"]["strelka"]["config"]),
        run_extra=""
    threads: 22
    wrapper:
        "0.65.0/bio/strelka/germline"

rule vcf_to_bcf:
    input:
        "{variants}.vcf.gz"
    output:
        "{variants}.output.bcf"
    log:
        "logs/bcftools/to-bcf/{variants}.log"
    params:
        "-O b -f PASS"
    wrapper:
        "0.60.0/bio/bcftools/view"


rule concat_somatic:
    input:
        calls=expand("results/strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf", type=["snvs","indels"]),
        indices=expand("results/strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf.csi", type=["snvs","indels"])
    output:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.bcf"
    log:
        "bcftools/concat-somatic/{sample}.log"
    params:
        "-O b -a"
    wrapper:
        "0.60.0/bio/bcftools/concat"

rule get_tumor_from_somatic:
    input:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.bcf"
    output:
        "results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf"
    log:
        "logs/bcftools/view-TUMOR/{sample}.log"
    params:
        "-O b -s TUMOR"
    wrapper:
        "0.60.0/bio/bcftools/view"

rule reheader_germline:
    input:
        vcf="{germline}/variants.output.bcf",
        samples="config/newsamples.txt"
    output:
        "{germline}/variants.reheader.bcf"
    log:
        "logs/bcftools/reheader/{germline}.log"
    params:
        extra="",
        view_extra="-O b"
    wrapper:
        "0.60.0/bio/bcftools/reheader"

rule concat_variants:
    input:
        # germline=get_germline_variants,
        # index_g=get_germline_variants_index,
        # somatic="results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf",
        # index_s="results/strelka/somatic/{sample}/results/variants/somatic.complete.tumor.bcf.csi",
        calls=lambda w: get_pair_variants(w, index=False),
        index=lambda w: get_pair_variants(w, index=True)
        #db="resources/snpEff/GRCh38.86"
    output:
        "results/strelka/merged/{sample}/all_variants.bcf"
    log:
        "bcftools/concat-all/{sample}.log"
    params:
        extra="-O b"
    wrapper:
        "0.64.0/bio/bcftools/concat"

rule annotate_strelka:
    input:
        calls="results/strelka/merged/{sample}/all_variants.bcf",
        db="resources/snpEff/GRCh38.86"
    output:
        calls="results/strelka/merged/{sample}/all_variants.annotated.bcf"
    log:
        "logs/snpeff/{sample}.log"
    params:
        #data_dir=lambda w, input: Path(input.db).parent.resolve(),
        #release=config["ref"]["snpeff"],
        extra="-Xmx4g -nodownload"
    wrapper:
        "0.64.0/bio/snpeff/annotate"

rule preprocess_variants:
    input:
        variants="{variants}.bcf"
    output:
        "{variants}.prepy.bcf"
    params:
        extra="-L --somatic",
        genome="resources/genome.fasta",
    log:
        "logs/prepy/{variants}.log"
    threads: 2
    wrapper:
        "0.60.0/bio/hap.py/pre.py"

# rule norm_vcf:
#     input:
#         "{prefix}.bcf",
#         genome="resources/genome.fasta"
#     output:
#         "{prefix}.norm.bcf"
#     log:
#         "logs/bcftools/norm/{prefix}.log"
#     params:
#         "-f {} -O b".format("resources/genome.fasta")  # optional parameters for bcftools norm (except -o)
#     wrapper:
#         "0.64.0/bio/bcftools/norm"
