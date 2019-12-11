rule arriba:
    input:
        bam="star/{sample}/Aligned.out.bam"
    output:
        fusions="arriba/{sample}.fusions.tsv",
        discarded="arriba/{sample}.fusions.discarded.tsv"
    params:
        genome=config["reference"]["genome"],
        annotation=config["reference"]["annotation"],
        blacklist=config["fusion"]["arriba"]["blacklist"]
    log:
        "logs/arriba/{sample}.log"
    conda:
        "../envs/arriba.yaml"
    threads:
        1
    shell:
        "arriba -x {input.bam} -O {output.discarded} -o {output.fusions} "
        "-a {params.genome} -g {params.annotation} -b {params.blacklist} -T -P > {log} 2>&1" #TODO: wrapper

## TODO: Update
#rule fusioncatcher:
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
