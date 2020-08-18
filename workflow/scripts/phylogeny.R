library(phangorn)

## read the variant matrix
mat = read.table(snakemake@input[["matrix"]], header=T)
## update colnames
colnames(mat) <-c("A","B","C","D","G","H","E")
## binarization
print(mat)
mat[mat=="."] <- 0
mat[mat>0] <- 1
print(mat)
## transpose matrix
mat_transposed <- t(mat)
## compute distance matrix
dist_mat <- dist.gene(mat_transposed)
## create nj-tree
njtree <- NJ(dist_mat)
## plot
#pdf(snakemake@output[["pdf"]])
pdf("test.pdf")
plot(njtree,"unrooted", main="NJ")
dev.off()
