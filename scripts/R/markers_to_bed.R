library(Seurat)
library(EnsDb.Mmusculus.v79)
library(rtracklayer)

markers <- commandArgs(trailingOnly = TRUE)[1]
bed.out <- commandArgs(trailingOnly = TRUE)[2]

commandArgs(trailingOnly = TRUE)[1]

markers     <- read.csv2(markers)
markers.pos <- markers[markers$p_val_adj < 0.05 & markers$avg_logFC > 0,]



genes.coords <- genes(EnsDb.Mmusculus.v79)
genes.coords <- genes.coords[genes.coords$gene_biotype %in% c("protein_coding","lincRNA")]


lapply(levels(markers$cluster),function(x){
  markers.x <- markers[markers$cluster == x,]
  markers.x <- head(markers.x,100)
  genes.x   <- genes.coords[genes.coords$symbol %in% markers.x$closest_gene,]
  rtracklayer::export(object = genes.x, con=paste0(bed.out,"/",x,"_marker_genes.bed"),format = 'bed')  
})
