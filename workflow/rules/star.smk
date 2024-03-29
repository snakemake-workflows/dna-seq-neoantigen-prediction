rule STAR_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf"
    output:
        directory("resources/STAR_index")
    params:
        sjdb_overhang="100",
        extra=""#"--sjdbGTFfile {} --sjdbOverhang 100".format(config["ref"]["annotation"])
    log:
        "logs/star/index.log"
    threads: 32
    cache: True
    wrapper:
        "0.42.0/bio/star/index"

rule align:
    input:
        "resources/STAR_index",
        fq1=lambda wc: units.loc[(wc.sample, "RNA"), "fq1"],
        fq2=lambda wc: units.loc[(wc.sample, "RNA"), "fq2"]
    output:
        # see STAR manual for additional output files
        "results/star/{sample}/Aligned.out.bam",
        "results/star/{sample}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="STAR_index",
        # optional parameters - designed to get chimeric alignments for fusion detection
        extra="--outSAMtype BAM Unsorted --chimSegmentMin 10 --chimOutType WithinBAM SoftClip "
            " --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 "
            " --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 "
            "--quantMode GeneCounts --sjdbGTFfile resources/genome.gtf {}".format(
            config["params"]["star"])
    threads: 8
    wrapper:
        "0.42.0/bio/star/align"
