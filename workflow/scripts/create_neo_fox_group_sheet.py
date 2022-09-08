import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

HLA_SUFFIXES_REGEX = r"[NLSCAQ]?"

# allowed loci according to NeoFox input data documentation:
# https://neofox.readthedocs.io/en/latest/03_01_input_data.html#file-with-patient-data
# * mhcIAlleles: comma separated MHC I alleles of the patient for HLA-A, HLA-B and
#   HLA-C. If homozygous, the allele should be added twice.
# * mhcIIAlleles: comma separated MHC II alleles of the patient for HLA-DRB1, HLA-DQA1,
#   HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.
ALLOWED_LOCI = {
    "A",
    "B",
    "C",
    "DRB1",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
}

mhc_alleles = pd.read_csv(
        snakemake.input.hla_la_bestguess,
        sep="\t",
    )
# the Allele column can contain multiple ";"-separated entries for the
# same locus -- NeoFox does a hard assertion that only two alleles per
# gene exist, so we chose to only ever keep the first of such multiple
# possibilities
mhc_alleles.loc[:, "Allele"] = mhc_alleles["Allele"].str.split(pat=";")
mhc_alleles = mhc_alleles.explode(["Allele"]).drop_duplicates(subset=["Locus", "Chromosome"])
mhc_alleles = mhc_alleles[mhc_alleles["Locus"].isin(ALLOWED_LOCI)]
mhc_alleles.loc[:, "Allele"] = (
    mhc_alleles["Allele"]
    .str
    .replace(
        r"([A-Z]+\d?)\*(\d+):(\d+)(:\d+)*G?(" + HLA_SUFFIXES_REGEX + r")",
        r"HLA-\1*\2:\3\5",
        regex=True,
    )
)
# the multiple ";"-separated entries from above can be identical after reducing
# to the allele group (1st number) and specific HLA protein (2nd number)
mhc_alleles = mhc_alleles.drop_duplicates(subset=["Chromosome", "Allele"])

mhc_one_alleles = ",".join( mhc_alleles.loc[ mhc_alleles["Locus"].str.len() == 1, "Allele"] )
mhc_two_alleles = ",".join( mhc_alleles.loc[ mhc_alleles["Locus"].str.len() > 1, "Allele"] )

patient_info = pd.DataFrame(
    data={
        "identifier": [ snakemake.params.group.name ],
        "tumorType": [ snakemake.params.group["tumorType"] ],
        "mhcIAlleles": [ mhc_one_alleles ],
        "mhcIIAlleles": [ mhc_two_alleles ],
    }
)

patient_info.to_csv(snakemake.output.group_sheet, sep="\t", quoting=3, index=False)
