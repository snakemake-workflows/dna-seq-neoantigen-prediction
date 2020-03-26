rule freebayes:
    input:
        ref=config["reference"]["genome"],
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

rule delly:
    input:
        ref=config["reference"]["genome"],
        samples=get_paired_bams,
        index=get_paired_bais,
    output:
        "results/candidate-calls/{pair}.delly.bcf"
    params:
        extra=config["params"].get("delly", "")
    log:
        "results/log/delly/{pair}.log"
    threads: 1 # Should number of threads scale by samples?
    wrapper:
        "0.42.0/bio/delly"
