rule kallisto_quant:
    input:
        fastq=get_quant_reads_input,
        index="resources/kallisto/transcripts.idx",
    output:
        directory("results/kallisto/{group}.{tumor_alias}"),
    params:
        extra=kallisto_params,
    log:
        "results/logs/kallisto/quant/{group}.{tumor_alias}.log",
    wrapper:
        "0.60.1/bio/kallisto/quant"


rule STAR_align:
    input:
        "resources/STAR_index",
        fq1=lambda wc: units.loc[(wc.sample, "RNA"), "fq1"],
        fq2=lambda wc: units.loc[(wc.sample, "RNA"), "fq2"],
        gtf="resources/genome.gtf",
    output:
        # see STAR manual for additional output files
        "results/star/{sample}/Aligned.out.bam",
        "results/star/{sample}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}.log",
    params:
        # path to STAR reference genome index
        index="STAR_index",
        # optional parameters - designed to get chimeric alignments for fusion detection
        extra=lambda wc, input: "--quantMode GeneCounts --sjdbGTFfile {} {}".format(
            input.gtf, config["params"]["star"]
        ),
    threads: 8
    wrapper:
        "0.42.0/bio/star/align"


rule arriba:
    input:
        bam="results/star/{sample}/Aligned.out.bam",
        genome="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        fusions="results/arriba/{sample}.fusions.tsv",
        discarded="results/arriba/{sample}.fusions.discarded.tsv",
    params:
        blacklist=config["fusion"]["arriba"]["blacklist"],
        extra=config["fusion"]["arriba"]["params"],
    log:
        "results/logs/arriba/{sample}.log",
    threads: 1
    wrapper:
        "0.60.1/bio/arriba"


## TODO: Update
# rule fusioncatcher:
#    input:
#        fq1=lambda w: units.loc[(w.sample, "RNA"), "fq1"],
#        fq2=lambda w: units.loc[(w.sample, "RNA"), "fq2"]
#    output:
#        directory("fusioncatcher/{sample}")
#    params:
#        extra="-T tmp -d ../../fusioncatcher_data"
#    log:
#        "logs/fusioncatcher/{sample}.log"
#    threads:
#        8
#    shell:
#        "fusioncatcher -i {input.fq1},{input.fq2} -o {output} {params.extra} -p {threads} > {log}"
