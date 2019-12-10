import pandas as pd
import glob

files = glob.glob("out/tables/*.filtered.tsv")

dfs, dfs_id = dict(), []

for f in files: 
	df = pd.read_csv(f, sep = '\t') 
	dfs[f] = df
	df = df[["ID_tumor"]]
	dfs_id.append(df)

ids = pd.concat(dfs_id)
ids=ids.groupby(ids.columns.tolist(), as_index = False).size().reset_index().rename(columns = {0:"occurrences"})
ids.to_csv("out/tables/candidate_occurrences.tsv", sep = '\t', index = False)

for k, v in dfs.items():
	df = v.merge(ids, on = "ID_tumor")
	df.to_csv(k.replace("filtered.tsv", "filtered.counts.csv"), sep = ',', index = False)
