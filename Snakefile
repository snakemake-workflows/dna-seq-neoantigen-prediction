import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")


normals = dict()
tumors = units[units.condition == "tumor"]["sample"].tolist()

for t in tumors:
    typ = units[units["sample"] == t]["type"].iloc[0]
    normals[t] = units[units.type == typ][units.condition == "healthy"]["sample"].iloc[0]

CHROMOSOMES = []
for i in range(1,22):
    CHROMOSOMES.append("chr"+str(i))
CHROMOSOMES.extend(["chrX"])

def allinput(wildcards):
    ret = []
    for t in tumors:
        ret.extend(expand(["microphaser/info/{tumor}_{normal}/filtered/{tumor}_{normal}.tsv", "microphaser/fasta/{tumor}_{normal}/filtered/{tumor}_{normal}.{chrom}.mt.fa", "microphaser/fasta/{tumor}_{normal}/filtered/{tumor}_{normal}.{chrom}.wt.fa"], tumor = t, normal = normals[t], chrom = CHROMOSOMES))
    return ret

rule all:
    input:
        expand(["bwa/{sample}.rmdup.bam","bwa/{sample}.rmdup.bam.bai"], sample = samples["sample"]),
        allinput



##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
#include: "rules/trim.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/microphaser.smk"
include: "rules/MHC_binding.smk"
