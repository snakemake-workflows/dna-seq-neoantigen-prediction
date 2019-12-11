rule strelka_somatic:
    input:
        normal=get_normal_bam,
        normal_index=get_normal_bai,
        tumor="bwa/{sample}.rmdup.bam",
        tumor_index="bwa/{sample}.rmdup.bam.bai"
    output:
        "strelka/somatic/{sample}/results/variants/somatic.snvs.vcf.gz",
        "strelka/somatic/{sample}/results/variants/somatic.indels.vcf.gz"
    log:
        "log/calling/strelka_somatic/{sample}.log"
    params:
        ref=config["reference"]["genome"],
        callRegions=config["reference"]["call_regions"],
        extra=config["params"]["strelka"]
    conda:
        "../envs/strelka.yaml"
    threads: 22
    shell:
        "configureStrelkaSomaticWorkflow.py --normalBam {input.normal} --tumorBam {input.tumor} "
        "--referenceFasta {params.ref} --runDir strelka/somatic/{wildcards.sample} "
        "--callRegions {params.callRegions} {params.extra} "
        "&& strelka/somatic/{wildcards.sample}/runWorkflow.py -m local -j {threads}"

rule strelka_germline:
    input:
        bam="bwa/{normal}.rmdup.bam",
        normal_index="bwa/{normal}.rmdup.bam.bai"
    output:
        "strelka/germline/{normal}/results/variants/variants.vcf.gz"
    log:
        "log/calling/strelka_germline/{normal}.log"
    params:
        ref=config["reference"]["genome"],
        callRegions=config["reference"]["call_regions"],
        extra=config["params"]["strelka"]
    conda:
        "../envs/strelka.yaml"
    threads: 22
    shell:
        "configureStrelkaGermlineWorkflow.py --bam {input.bam} "
        "--referenceFasta {params.ref} --runDir strelka/germline/{wildcards.normal} "
        "--callRegions {params.callRegions} {params.extra} "
        "&& strelka/germline/{wildcards.normal}/runWorkflow.py -m local -j {threads}"

rule vcf_to_bcf:
    input:
        "{variants}.vcf.gz"
    output:
        "{variants}.output.bcf"
    params:
        "-O u -f PASS"
    wrapper:
        "0.31.1/bio/bcftools/view"

rule compress_bcf:
    input:
        "{variants}.bcf"
    output:
        "{variants}.bcf.gz"
    params:
        "-O b -s TUMOR"
    wrapper:
        "0.31.1/bio/bcftools/view"

rule index_bcf:
    input:
        "{variants}.bcf"
    output:
        "{variants}.bcf.csi"
    params:
        ""
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools index {input}"


rule concat_somatic:
    input:
        calls=expand("strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf", type=["snvs","indels"]),
        indices=expand("strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf.csi", type=["snvs","indels"])
    output:
        "strelka/somatic/{sample}/results/variants/somatic.complete.bcf"
    params:
        "-O b -a"
    wrapper:
        "0.31.1/bio/bcftools/concat"

rule clean_germline:
    input:
        bcf="{germline}/variants.output.bcf"
    output:
        "{germline}/variants.reheader.bcf"
    params:
        extra=""
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools reheader {params.extra} -s ref/newsamples.txt {input.bcf} > {output}"

rule concat_variants:
    input:
        germline=get_germline_variants,
        index_g=get_germline_variants_index,
        somatic="strelka/somatic/{sample}/results/variants/somatic.complete.bcf",
        index_s="strelka/somatic/{sample}/results/variants/somatic.complete.bcf.csi"
    output:
        "strelka/merged/{sample}/all_variants.bcf"
    params:
        assembly=config["reference"]["assembly"],
        extra="-O v"
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools concat {params.extra} -a {input.somatic} {input.germline} | snpEff -Xmx4g {params.assembly} - | bcftools view -O u - > {output}"
