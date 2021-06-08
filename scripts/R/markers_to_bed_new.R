library(argparse)

parser <- ArgumentParser()

parser$add_argument("-i", "--input", type="character", default='', 
                    help="path to csv2 with markers")

parser$add_argument("-n", "--nmarkers",type="integer", default=100, 
                    help="number of markers to export")

parser$add_argument("-o", "--out", type="character", default="foo", 
                    help="folder to output the bed files to")

args      <- parser$parse_args()

library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ensembldb)
library(rtracklayer)


markers         <- read.csv2(file = args$input)
markers$cluster <- gsub(" ","_",markers$cluster)
markers.pos     <- markers[markers$p_val < 0.05 & markers$avg_logFC > 0,]

conversion.table <- as.data.frame(genes(EnsDb.Mmusculus.v79))[,c("gene_name","symbol")]

genes.promoters        <- promoters(EnsDb.Mmusculus.v79)
genes.promoters        <- genes.promoters[genes.promoters$tx_biotype %in% c("protein_coding","lincRNA")]
genes.promoters$symbol <- conversion.table[genes.promoters$gene_id,"symbol",drop=TRUE]

genes.genes <- genes(EnsDb.Mmusculus.v79)
genes.genes <- genes.genes[genes.genes$gene_biotype %in% c("protein_coding","lincRNA"),]


clusters_to_use <- unique(markers$cluster)


print(getwd())
dir.create(path = args$out, recursive = TRUE,showWarnings = TRUE)


lapply(clusters_to_use,function(x){
  markers.x <- markers.pos[markers.pos$cluster == x,]
  markers.x <- head(markers.x,args$nmarkers)
  promoters.x <- genes.promoters[genes.promoters$symbol %in% markers.x$closest_gene | genes.promoters$symbol %in% markers.x$gene ]
  rtracklayer::export(object = promoters.x,con = paste0(args$out,"/",x,"_promoters.bed"))
})


lapply(clusters_to_use,function(x){
  markers.x <- markers.pos[markers.pos$cluster == x,]
  markers.x <- head(markers.x,args$nmarkers)
  genes.x <- genes.genes[genes.genes$symbol %in% markers.x$closest_gene | genes.genes$symbol %in% markers.x$gene]
  rtracklayer::export(object = genes.x,con = paste0(args$out,"/",x,"_genes.bed"))
})




