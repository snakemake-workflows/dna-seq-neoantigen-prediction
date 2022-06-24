import pandas as pd

# to get alleles that netMHCpan can handle, use its -listMHC option
hlaI = ["A","B","C", "E", "G"]

# to get alleles that netMHCIIpan can handle, use its -list option
hlaII = ["DRB1", "DRB3", "DRB4", "DRB5", "DPA1", "DPB1", "DQA1", "DQB1"]

hlas = pd.read_csv(snakemake.input[0], sep='\t')

hlasI = hlas[hlas.Locus.isin(hlaI)]
hlasI["Allele"]="HLA-" + hlasI.Allele.str.split(":", expand=True)[[0,1]].apply(lambda x: ''.join(x), axis=1).str.replace('*','')
hlasI = hlasI[["Allele"]].drop_duplicates()
hlasI.to_csv(snakemake.output[0], sep='\t', index=False)

hlasII = hlas[hlas.Locus.isin(hlaII)]
hlasII["HLA"] = hlasII.Locus.str[0:2]
hlasII["Allele"] = hlasII.Allele.str.split(":", expand=True)[[0,1]].apply(lambda x: ''.join(x), axis=1).str.replace('*','')

hlasII = pd.DataFrame("HLA-" + hlasII.groupby(["HLA","Chromosome"])["Allele"].apply(lambda x: "-".join(x)).reset_index()["Allele"]).drop_duplicates()
hlasII.loc[hlasII.Allele.str.contains("DRB"), "Allele"] = hlasII[hlasII.Allele.str.contains("DRB")]["Allele"].str.replace("HLA-DRB1","DRB1_")
hlasII.to_csv(snakemake.output[1], sep='\t', index=False)
