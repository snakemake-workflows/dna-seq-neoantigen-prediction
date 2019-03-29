def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"])
