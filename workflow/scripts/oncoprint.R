# log to file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ComplexHeatmap)
library(ggplot2)

table = read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mat = as.matrix(table)
mat = t(mat)
## remove "full" lines
mat = mat[rowSums(mat == "") > 0,]
## remove single lines
mat = mat[rowSums(mat != "") > 0,]

col = c(SNV = "blue", INDEL = "red")

alter_fun = list(
        SNV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
            gp = gpar(fill = col["SNV"], col = NA)),
        INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
            gp = gpar(fill = col["INDEL"], col = NA))
    )


heatmap_legend_param = list(title = "Alterations", at = c("SNV", "INDEL"), 
        labels = c("SNV", "INDEL"))
if (ncol(mat) > 1 ){
    mat <- mat[order(apply(mat, 1, function(row) sum(row != "")), decreasing = T), ]
}

i = 0
c = 0
matlist = list()
while (i + 2000 < nrow(mat)) {
    m <- mat[(i + 1):(i + 2000), , drop=FALSE]
    rows_matrix <- nrow(m)
    height_plot <- (rows_matrix/5)
    if (height_plot < 4) {
        height_plot <- 4
    }
    pdf(file = sub("0.pdf", paste(c, ".pdf", sep=''), snakemake@output[[1]]), height=height_plot)
    if (rows_matrix > 0) {
        oncoprint <- oncoPrint(m,
            alter_fun = alter_fun, col = col, 
            remove_empty_columns = FALSE, remove_empty_rows = TRUE,
            pct_side = "right", row_names_side = "left",
            show_column_names=T,
            column_title = "OncoPrint", heatmap_legend_param = heatmap_legend_param)
        draw(oncoprint, newpage=F)
    }
    dev.off()
    i = i + 2000
    c = c + 1
}