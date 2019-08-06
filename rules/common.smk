def is_single_end(sample, typ):
    return pd.isnull(units.loc[(sample, typ), "fq2"])

def get(sample, typ):
    return [units.loc[(sample, typ), "fq1"], units.loc[(sample, typ), "fq2"]]

def get_seperate(sample, typ, group):
    return units.loc[(sample, typ), "fq" + str(group)]
