import sys

sys.stderr = open(snakemake.log[0], "w")

import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pysam import VariantFile

variant_df = pd.read_csv(snakemake.input[0], sep="\t").fillna(0.0)
variant_df = variant_df[
    ["CHROM", "POS"] + [c for c in variant_df.columns if c.endswith("Freq")]
]
## tidy data - for facet plot
tidy_df = variant_df.melt(id_vars=["CHROM", "POS"], var_name="Sample", value_name="VAF")
g = sns.FacetGrid(tidy_df, col="Sample")
g = g.map(sns.distplot, "VAF", hist=False)
g.savefig(snakemake.output["facet"])
plt.close()
## pairplot
sns.pairplot(variant_df.drop(["CHROM", "POS"], axis=1), diag_kind="kde")
plt.savefig(snakemake.output["pairplot"])
plt.close()


def overlap_pct(x, y, **kws):
    n = 0
    for i in range(0, len(x)):
        if (x[i] > 0) & (y[i] > 0):
            n += 1
    overlap = n / len([e for e in x if e > 0])
    ax = plt.gca()
    ax.annotate("Shared Fraction: {:.2f}".format(overlap), xy=(0.2, 0.4))
    ax.annotate("Shared Variants: {}".format(n), xy=(0.2, 0.6))


def variants(x, **kws):
    positive = len([e for e in x if e > 0])
    ax = plt.gca()
    ax.annotate("#Variants: {}".format(positive), xy=(0.2, 0.5), xycoords=ax.transAxes)


g = sns.PairGrid(variant_df.drop(["CHROM", "POS"], axis=1))

g.map_offdiag(overlap_pct)
g.map_diag(variants)
g.savefig(snakemake.output["grid"])
plt.close()

# def neg_overlap_pct(x, y, **kws):
#    n = 0
#    for i in range(0,len(x)):
#        if (x[i] == 0) & (y[i] == 0):
#            n += 1
#    overlap=n/len([e for e in x if e == 0])
#    ax = plt.gca()
#    ax.annotate("Shared Fraction: {:.2f}".format(overlap), xy=(.2, .4))
#    ax.annotate("Shared missing variants: {}".format(n), xy=(.2, .6))

# def neg_variants(x, **kws):
#    zero = len([e for e in x if e == 0])
#    ax = plt.gca()
#    ax.annotate("#Missing Variants: {}".format(zero), xy=(.2, .5),
# xycoords=ax.transAxes)

# g = sns.PairGrid(variant_df.drop(["CHROM", "POS"], axis=1))

# g.map_offdiag(neg_overlap_pct)
# g.map_diag(neg_variants)
# g.savefig("plots/Mssing_Variant_table.pdf")
# plt.close()

for c in variant_df.drop(["CHROM", "POS"], axis=1).columns:
    sns.distplot(variant_df[[c]][variant_df[c] > 0])
    plt.title(c)
    plt.savefig("plots/positive_" + c + ".distplot.pdf")
    plt.close()
    sns.distplot(variant_df[[c]])
    plt.savefig("plots/all_" + c + ".distplot.pdf")
    plt.close()
