def get_normal_bam(wildcards):
    return(expand("bwa/{normal}.rmdup.bam", normal=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_normal_bai(wildcards):
    return(expand("bwa/{normal}.rmdup.bam.bai", normal=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

rule strelka_somatic:
    input:
        normal=get_normal_bam,#"bwa/{normal}.rmdup.bam",
        normal_index=get_normal_bai,#"bwa/{normal}.rmdup.bam.bai",
        tumor="bwa/{sample}.rmdup.bam",
        tumor_index="bwa/{sample}.rmdup.bam.bai"
    output:
#        directory("strelka/somatic/{sample}")
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
#        directory("strelka/germline/sample")
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
        "{variants}.bcf.gz"
    output:
        "{variants}.bcf.gz.csi"
    params:
        ""
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools index {input}"


rule concat_somatic:
    input:
        calls=expand("strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf.gz", type=["snvs","indels"]),
        indices=expand("strelka/somatic/{{sample}}/results/variants/somatic.{type}.output.bcf.gz.csi", type=["snvs","indels"])
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

def get_germline(wildcards):
    return(expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

def get_germline_index(wildcards):
    return(expand("strelka/germline/{germline}/results/variants/variants.reheader.bcf.gz.csi", germline=samples[samples["sample"] == wildcards.sample]["matched_normal"]))

rule concat_variants:
    input:
        germline=get_germline,#"strelka/germline/{sample}/results/variants/variants.reheader.bcf.gz",
        index_g=get_germline_index,#"strelka/germline/{sample}/results/variants/variants.reheader.bcf.gz.csi",
        somatic="strelka/somatic/{sample}/results/variants/somatic.complete.bcf.gz",
        index_s="strelka/somatic/{sample}/results/variants/somatic.complete.bcf.gz.csi"
    output:
        "strelka/merged/{sample}/all_variants.bcf"
    params:
        assembly=config["reference"]["assembly"],
        extra="-O v"
    conda:
        "../envs/variant_handling.yaml"
    shell:
        "bcftools concat {params.extra} -a {input.somatic} {input.germline} | snpEff -Xmx4g {params.assembly} - | bcftools view -O u - > {output}"
