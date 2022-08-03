import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

# To know which alleles netMHCpan can handle, use its -listMHC option.
HLAI = {"A", "B", "C", "E", "G"}

# To know which alleles netMHCIIpan can handle, use its -list option.
# DRB alleles need to be formatted differently from DP and DQ alleles,
# so we specify them separately.
DRB = {"DRB1", "DRB3", "DRB4", "DRB5"}
ALPHA_BETA = {"DPA1", "DPB1", "DQA1", "DQB1"}

hlas = pd.read_csv(snakemake.input.hla_la_bestguess, sep="\t")

# the Allele column can contain multiple ";"-separated entries for the
# same locus
hlas.loc[:, "Allele"] = hlas["Allele"].str.split(pat=";")
hlas["alternative"] = hlas["Allele"].apply(lambda x: range(len(x)))
hlas = hlas.explode(["Allele", "alternative"])

# reformat to netMHCpan allele list format:
# * https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list
# it needs to be in the format of the first column of the above list, as explained in
# the "Instructions" tab under "MHC SELECTION" point "2)" at:
# * https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
hlaI_alleles = (
    hlas.loc[hlas["Locus"].isin(HLAI), "Allele"]
    .str.replace(r"([A-Z])\*(\d+):(\d+)(:\d+)*G?(N?)", r"HLA-\1\2:\3\5", regex=True)
    .drop_duplicates()
)
hlaI_alleles.to_csv(snakemake.output.hlaI, sep="\t", index=False, header=False)

# reformat to netMHCIIpan allele list format:
# * https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.1/alleles_name.list
# contrary to the format in that list, alleles actually need to be formatted like this,
# with <gene>s found in the HLA-LA "Locus" column and syntax for the sub-numbering (only
# the 1st and 2nd sub-number are used) according to the official nomenclature (see:
# https://ars.els-cdn.com/content/image/1-s2.0-S0006497120405555-gr2.jpg ):
# * DRB alleles: "<gene>_<allele_group><specific_HLA_protein>"
# * DP and DQ alleles (alpha means A and beta means B in the gene name, for example DPA):
#   "HLA-<alpha_gene><allele_group><specific_HLA_protein>-<beta_gene><allele_group><specific_HLA_protein>"
# This format was determined by manually selecting combinations above the
# "type a list of molecules names" field of the "Submission" tab at:
# * https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.1

# TODO: check whether Jan's previous parsing of DRB alleles into this format is necessary:
# * example: DRB1_1501-DRB30101-DRB40301N
# * "DRB1_<allele_group><specific_HLA_protein>-DRB3<allele_group><specific_HLA_protein>-DRB4<allele_group><specific_HLA_protein>"
drb_alleles = hlas.loc[hlas["Locus"].isin(DRB)]
hlaII_alleles = (
    drb_alleles["Allele"]
    .str.replace(r"([A-Z]+\d)\*(\d+):(\d+)(:\d+)*G?(N?)", r"\1_\2\3\5", regex=True)
    .drop_duplicates()
)

# handle alleles where a combination of alpha and beta always exists
alpha_beta_alleles = hlas.loc[hlas["Locus"].isin(ALPHA_BETA)]
alpha_beta_alleles.loc[:, "Allele"] = alpha_beta_alleles.Allele.str.replace(
    r"([A-Z]+\d)\*(\d+):(\d+)(:\d+)*G?(N?)", r"\1\2\3\5", regex=True
)
# we need a variable to group alpha and beta of the same gene combination together
alpha_beta_alleles["gene_group"] = alpha_beta_alleles["Locus"].str[0:2]
# we need to handle cases where we had multiple allele entries in an
# alpha or beta locus, adding in a duplicate of the corresponding locus
select_mult_all = alpha_beta_alleles["alternative"] > 0
select_dpq_loci = alpha_beta_alleles["Locus"].str.startswith("D[PQ]")
mult_all_per_loc_selection = select_mult_all & select_dpq_loci
alleles_to_duplicate = alpha_beta_alleles.loc[
    mult_all_per_loc_selection,
    ["Locus", "Chromosome", "alternative"],
].replace(regex={"Locus": {"(D[PQ])A(\d+)": r"\1B\2", "(D[PQ])B(\d+)": r"\1A\2"}})
alleles_to_insert = alleles_to_duplicate.merge(
    alpha_beta_alleles.drop("alternative", axis="columns"),
    on=["Locus", "Chromosome"],
    how="left",
)
alpha_beta_alleles = pd.concat(
    [alpha_beta_alleles, alleles_to_insert]
).drop_duplicates()

hlaII_alleles = hlaII_alleles.append(
    alpha_beta_alleles.groupby(["gene_group", "Chromosome", "alternative"])["Allele"]
    .transform(lambda x: f"HLA-{'-'.join(x)}")
    .drop_duplicates()
)

hlaII_alleles.to_csv(snakemake.output.hlaII, sep="\t", index=False, header=False)
