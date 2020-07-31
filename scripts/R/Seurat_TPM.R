library(Seurat)
library(Matrix.utils)
library(EnsDb.Mmusculus.v79)

args <- commandArgs(trailingOnly=TRUE)

input.seurat <- args[1]
out.csv      <- args[2]

# Load and aggregate seuray
brain <- readRDS(input.seurat)
brain.agg <- aggregate.Matrix(x= t(brain[['RNA']]@counts), groupings = brain@active.ident, FUN='sum')

genes <- ensembldb::genes(EnsDb.Mmusculus.v79)

# Load annotations - lengths
genes.lengths           <- data.frame(length = width(genes), gene = genes$symbol,ensg=genes$gene_id)
rownames(genes.lengths) <- genes.lengths$ensg

genes.lengths.agg <- aggregate(genes.lengths$length, by= list(gene=genes.lengths$ensg),FUN='mean')
genes.lengths.agg$symbol <- genes.lengths$gene[genes.lengths.agg$gene]

# Rename the matrix
common.genes <- intersect(genes.lengths$gene,colnames(brain.agg))

# Filter common genes
genes.lengths <- genes.lengths[genes.lengths$gene %in% common.genes,]
brain.agg     <- brain.agg[,colnames(brain.agg)%in% common.genes]


# Change gene symbol to ENSG in matrix
rownames(genes.lengths.agg) <- make.unique(as.character(genes.lengths.agg$symbol))
colnames(brain.agg)         <- as.character(genes.lengths.agg[colnames(brain.agg),'gene'])
rownames(genes.lengths.agg) <- as.character(genes.lengths.agg$gene)

brain.agg.RPK <-  apply(brain.agg,1,function(x){ x / (genes.lengths.agg[names(x),'x']/1000)  })
brain.agg.TPM <-  as.data.frame(apply(brain.agg.RPK,2,function(x){x/(sum(x)/1e6)}))

write.table(brain.agg.TPM,out.csv,quote=FALSE)