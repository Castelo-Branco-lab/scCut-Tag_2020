library(Seurat)
library(Matrix.utils)
library(EnsDb.Mmusculus.v79)

args <- commandArgs(trailingOnly=TRUE)

input.seurat <- args[1]
out.csv      <- args[2]

brain <- readRDS(input.seurat)
brain.agg <- aggregate.Matrix(x= t(brain[['RNA']]@counts), groupings = brain@active.ident, FUN='sum')

genes <- ensembldb::genes(EnsDb.Mmusculus.v79)

genes.lengths <- data.frame(length = width(genes), gene = genes$symbol)
genes.lengths <- aggregate(genes.lengths$length, by= list(gene=genes.lengths$gene),FUN='mean')

common.genes <- intersect(genes.lengths$gene,colnames(brain.agg))

genes.lengths <- genes.lengths[genes.lengths$gene %in% common.genes,]
brain.agg     <- brain.agg[,colnames(brain.agg)%in% common.genes]

rownames(genes.lengths) <- genes.lengths$gene


brain.agg.RPK <-  apply(brain.agg,1,function(x){ x / (genes.lengths[names(x),'x']/1000)  })
brain.agg.TPM <- apply(brain.agg.RPK,2,function(x){x/(sum(x)/1e6)})

write.table(brain.agg.TPM,out.csv,quote=FALSE)