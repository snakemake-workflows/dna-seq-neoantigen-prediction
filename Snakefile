import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t').set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep='\t').set_index(["sample", "type"], drop=False)
#validate(units, schema="schemas/units.schema.yaml")



SAMPLES = samples["sample"]

CHROMOSOMES = []
for i in range(1,23):
    CHROMOSOMES.append("chr"+str(i))
CHROMOSOMES.extend(["chrX"])

def allinput(wildcards):
    ret = []
    tumors = samples[samples.condition == "tumor"].index.to_list()
    ret.extend(expand("results/{mhc}/{sample}.WES.tsv", sample=tumors, mhc=["netMHCpan", "netMHC2pan"]))
    return ret

rule all:
    input:
        expand("results/{mhc}/{sample}.WES.tsv", sample=samples[samples.condition == "tumor"]["sample"], mhc=["netMHCpan", "netMHC2pan"])



##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
#include: "rules/quant.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/microphaser.smk"
include: "rules/MHC_binding.smk"
