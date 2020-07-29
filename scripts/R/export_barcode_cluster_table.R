library(Seurat)
library(Rsamtools)

args = commandArgs(trailingOnly=TRUE)

brain     <- readRDS(args[1])
outfile   <- args[2]

cluster.annotations           <- data.frame(cluster=brain@active.ident,barcode=names(brain@active.ident))
colnames(cluster.annotations) <- paste0("#",colnames(cluster.annotations))



write.csv(x = cluster.annotations, file=args[2],row.names = FALSE,quote = FALSE,)
