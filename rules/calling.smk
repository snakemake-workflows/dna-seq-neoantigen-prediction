rule Strelka_somatic:
    input:
        normal="bwa/{normal}.rmdup.bam",
        index_n="bwa/{normal}.rmdup.bam.bai",
        tumor="bwa/{tumor}.rmdup.bam",
        index_t="bwa/{tumor}.rmdup.bam.bai"
    output:
        "strelka/{tumor}-{normal}/results/variants/somatic.snvs.vcf.gz",
        "strelka/{tumor}-{normal}/results/variants/somatic.indels.vcf.gz"
    log:
        "log/calling/strelka_somatic/{tumor}_{normal}.log"
    params:
        ref=config["reference"]["genome"],
        callRegions=config["reference"]["call_regions"],
        extra=config["params"]["strelka"]
    conda:
        "../envs/strelka.yaml"
    threads: 22
    shell:
        "configureStrelkaSomaticWorkflow.py --normalBam {input.normal} --tumorBam {input.tumor} "
        "--referenceFasta {params.ref} --runDir strelka/{wildcards.tumor}-{wildcards.normal} "
        "--callRegions {params.callRegions} {params.extra} "
        "&& strelka/{wildcards.tumor}-{wildcards.normal}/runWorkflow.py -m local -j {threads}"

rule strelka_germline:
    input:
        bam="bwa/{normal}.rmdup.bam",
        index_n="bwa/{normal}.rmdup.bam.bai"
    output:
        "strelka/{normal}/results/variants/variants.vcf.gz"
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
        "--referenceFasta {params.ref} --runDir strelka/{wildcards.normal} "
        "--callRegions {params.callRegions} {params.extra} "
        "&& strelka/{wildcards.normal}/runWorkflow.py -m local -j {threads}"

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
        calls=expand("strelka/{{tumor}}-{{normal}}/results/variants/somatic.{type}.output.bcf.gz", type=["snvs","indels"]),
        indices=expand("strelka/{{tumor}}-{{normal}}/results/variants/somatic.{type}.output.bcf.gz.csi", type=["snvs","indels"])
    output:
        "strelka/{tumor}-{normal}/results/variants/somatic.complete.bcf"
    params:
        "-O b -a"
    wrapper:
        "0.31.1/bio/bcftools/concat"

#rule clean_germline:
#    input:
#        bcf="strelka/{normal}/results/variants/variants.output.bcf"
#    output:
#        "strelka/{normal}/results/variants/variants.reheader.bcf"
#    params:
#        ""
#    conda:
#        "../envs/variant_handling.yaml"
#    shell:
#        "bcftools reheader {params} -s ref/newsamples.txt {input.bcf} > {output}"


#rule concat_variants:
#    input:
#        germline="strelka/{normal}/results/variants/variants.reheader.bcf.gz",
#        index_g="strelka/{normal}/results/variants/variants.reheader.bcf.gz.csi",
#        somatic="strelka/{tumor}-{normal}/results/variants/somatic.complete.bcf.gz",
#        index_s="strelka/{tumor}-{normal}/results/variants/somatic.complete.bcf.gz.csi"
#    output:
#        "strelka/{tumor}-{normal}/results/variants/all_variants.bcf"
#    params:
#        assembly=config["reference"]["assembly"],
#        extra="-Ov"
#    conda:
#        "../envs/variant_handling.yaml"
#    shell:
#        "bcftools concat {params.extra} -a {input.somatic} {input.germline} | snpEff -Xmx4g {params.assembly} - | bcftools view -Ou - > {output}"

rule preprocess_variants:
    input:
        "{variants}.bcf"
    output:
        "{variants}.prepy.bcf"
    params:
        extra="-L --somatic",
        ref=config["reference"]["genome"],
    conda:
        "../envs/prepy.yaml"
    shell:
        "pre.py {params.extra} -r {params.ref} {input} {output}"

