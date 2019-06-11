def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"])

def get(sample):
    return [units.loc[(sample), "fq1"], units.loc[(sample), "fq2"]]

def get_seperate(sample, group):
    return units.loc[(sample), "fq" + str(group)]
