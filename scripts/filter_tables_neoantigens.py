import pandas as pd

data = pd.read_csv(snakemake.input["table"], sep = '\t')

data = data[["ID_tumor","Gene_ID","Gene_Symbol","Chromosome","offset","freq_tumor",
"AminoAcid_Change","Peptide_tumor","NB_tumor","Rank_min_tumor","Aff_min_tumor",
"Top_ranked_HLA_tumor","Top_affinity_HLA_tumor","Peptide_normal","NB_normal",
"Rank_min_normal","Aff_min_normal","Top_ranked_HLA_normal","Top_affinity_HLA_normal"]]

data.columns = ["ID_tumor","Gene_ID","Gene_Symbol","Chromosome","Position","Frequency_tumor",
"AminoAcid_Change","Peptide_tumor","BindingHLAs_tumor","Rank_min_tumor","Affinity_min_tumor",
"Top_rank_HLA_tumor","Top_affinity_HLA_tumor","Peptide_normal","BindingHLAs_normal",
"Rank_min_normal","Aff_min_normal","Top_rank_HLA_normal","Top_affinity_HLA_normal"]

data = data[data.BindingHLAs_tumor > 0]
# data = data[(data.NB_normal.isna()) | (data.NB_normal == 0)]
data = data[(data.BindingHLAs_normal == 0)]

### Delete Stop-Codon including peptides
data = data[data.Peptide_tumor.str.count("x") == 0]
data = data[data.Peptide_tumor.str.count("X") == 0]

data.to_csv(snakemake.output[0], sep = '\t', index = False)
