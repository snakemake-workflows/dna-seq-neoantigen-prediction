rule freebayes:
    input:
        ref="resources/genome.fasta",
        # you can have a list of samples here
        samples=get_paired_bams
    output:
        "results/candidate-calls/{pair}.freebayes.bcf"
    log:
        "results/log/freebayes/{pair}.log"
    params:
        extra=config["params"].get("freebayes", ""),
        chunksize=100000
    threads: 60
    wrapper:
        "0.42.0/bio/freebayes"
