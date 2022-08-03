import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

# # read in available alleles


def read_allele_list(filename: str):
    with open(filename, "r") as alleles_in:
        alleles = set()
        for line in alleles_in:
            if not (line.startswith("#") or line == "\n"):
                alleles.add(line.strip())
        return alleles


HLA_SUFFIXES_REGEX = r"[NLSCAQ]?"

# netMHCpan alleles and loci
HLA_ONE_NET_MHC_ALLELES = read_allele_list(snakemake.input.mhc_one_alleles)
hla_one_net_mhc_alleles = pd.Series(list(HLA_ONE_NET_MHC_ALLELES))
HLA_ONE_LOCI = set(
    hla_one_net_mhc_alleles[hla_one_net_mhc_alleles.str.startswith("HLA-")]
    .str.replace(r"HLA-([A-Z])[\d:]+" + HLA_SUFFIXES_REGEX, r"\1", regex=True)
    .drop_duplicates()
)


# netMHCIIpan alleles and loci
HLA_TWO_NET_MHC_ALLELES = read_allele_list(snakemake.input.mhc_two_alleles)
hla_two_net_mhc_alleles = pd.Series(list(HLA_TWO_NET_MHC_ALLELES))
# DRB alleles need to be formatted differently from DP and DQ alleles,
# so we extract them separately.
DRB_LOCI = set(
    hla_two_net_mhc_alleles[hla_two_net_mhc_alleles.str.startswith("DRB")]
    .str.replace(r"(DRB\d)_\d+" + HLA_SUFFIXES_REGEX, r"\1", regex=True)
    .drop_duplicates()
)

ALPHA_BETA_LOCI = set(
    hla_two_net_mhc_alleles[hla_two_net_mhc_alleles.str.startswith("HLA-")]
    .str.replace(
        r"HLA-(D[A-Z]A\d)\d+"
        + HLA_SUFFIXES_REGEX
        + r"-(D[A-Z]B\d)\d+"
        + HLA_SUFFIXES_REGEX,
        r"\1_\2",
        regex=True,
    )
    .str.split("_")
    .explode()
    .drop_duplicates()
)

# read in alleles as determined by HLA-LA
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
hla_one_alleles = (
    hlas.loc[hlas["Locus"].isin(HLA_ONE_LOCI), "Allele"]
    .str.replace(
        r"([A-Z])\*(\d+):(\d+)(:\d+)*G?(" + HLA_SUFFIXES_REGEX + r")",
        r"HLA-\1\2:\3\5",
        regex=True,
    )
    .drop_duplicates()
)
hla_one_alleles[hla_one_alleles.isin(HLA_ONE_NET_MHC_ALLELES)].to_csv(
    snakemake.output.hlaI, sep="\t", index=False, header=False
)

# reformat to netMHCIIpan allele list format:
# * https://services.healthtech.dtu.dk/services/netMHCIIpan-4.1/alleles_name.list
# contrary to the format in that list, alleles actually need to be formatted like this,
# with <gene>s found in the HLA-LA "Locus" column and syntax for the sub-numbering (only
# the 1st and 2nd sub-number are used) according to the official nomenclature (see:
# http://www.hla.alleles.org/nomenclature/naming.html ):
# * DRB alleles: "<gene>_<allele_group><specific_HLA_protein>"
# * DP and DQ alleles (alpha means A and beta means B in the gene name, for example DPA):
#   "HLA-<alpha_gene><allele_group><specific_HLA_protein>-<beta_gene><allele_group><specific_HLA_protein>"
# This format was determined by manually selecting combinations above the
# "type a list of molecules names" field of the "Submission" tab at:
# * https://services.healthtech.dtu.dk/service.php?netMHCIIpan-4.1

# TODO: check whether Jan's previous parsing of DRB alleles into this format is necessary:
# * example: DRB1_1501-DRB30101-DRB40301N
# * "DRB1_<allele_group><specific_HLA_protein>-DRB3<allele_group><specific_HLA_protein>-DRB4<allele_group><specific_HLA_protein>"
drb_alleles = hlas.loc[hlas["Locus"].isin(DRB_LOCI)]
hla_two_alleles = (
    drb_alleles["Allele"]
    .str.replace(
        r"([A-Z]+\d)\*(\d+):(\d+)(:\d+)*G?(" + HLA_SUFFIXES_REGEX + r")",
        r"\1_\2\3\5",
        regex=True,
    )
    .drop_duplicates()
)

# handle alleles where a combination of alpha and beta always exists
alpha_beta_alleles = hlas.loc[hlas["Locus"].isin(ALPHA_BETA_LOCI)]
alpha_beta_alleles.loc[:, "Allele"] = alpha_beta_alleles.Allele.str.replace(
    r"([A-Z]+\d)\*(\d+):(\d+)(:\d+)*G?(" + HLA_SUFFIXES_REGEX + r")",
    r"\1\2\3\5",
    regex=True,
)
# we need a variable to group alpha and beta of the same gene combination together
alpha_beta_alleles.loc[:, "gene_group"] = alpha_beta_alleles["Locus"].str.replace(
    r"(D[A-Z])[AB]\d", r"\1", regex=True
)
# we need to handle cases where we had multiple allele entries in an
# alpha or beta locus, adding in a duplicate of the corresponding locus
alleles_to_duplicate = alpha_beta_alleles.loc[
    alpha_beta_alleles["alternative"] > 0,
    ["Locus", "Chromosome", "alternative"],
].replace(regex={"Locus": {"(D[A-Z])A(\d)": r"\1B\2", "(D[A-Z])B(\d)": r"\1A\2"}})
alleles_to_insert = alleles_to_duplicate.merge(
    alpha_beta_alleles.drop("alternative", axis="columns"),
    on=["Locus", "Chromosome"],
    how="left",
)
alpha_beta_alleles = pd.concat(
    [alpha_beta_alleles, alleles_to_insert]
).drop_duplicates()

hla_two_alleles = hla_two_alleles.append(
    alpha_beta_alleles.groupby(["gene_group", "Chromosome", "alternative"])["Allele"]
    .transform(lambda x: f"HLA-{'-'.join(x)}")
    .drop_duplicates()
)

hla_two_alleles[hla_two_alleles.isin(HLA_TWO_NET_MHC_ALLELES)].to_csv(
    snakemake.output.hlaII, sep="\t", index=False, header=False
)
