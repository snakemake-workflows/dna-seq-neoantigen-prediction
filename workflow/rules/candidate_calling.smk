rule freebayes:
    input:
        ref=config["reference"]["genome"],
        # you can have a list of samples here
        samples=get_paired_bams
    output:
        "candidate-calls/{group}.freebayes.bcf"
    log:
        "logs/freebayes/{group}.log"
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
        "candidate-calls/{group}.delly.bcf"
    params:
        extra=config["params"].get("delly", "")
    log:
        "logs/delly/{group}.log"
    threads: 1 # Should number of threads scale by samples?
    wrapper:
        "0.42.0/bio/delly"
