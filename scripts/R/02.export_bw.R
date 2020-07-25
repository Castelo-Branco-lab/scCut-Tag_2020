library(Seurat)
library(Signac)
require(GenomicRanges)
require(tools)

# Debug
#args <- c("clustering_snakemake/01.clustering/Seurat_5000_UMI.Rds","monocle")
  args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2){
  stop("Must provide 2 arguments: 1. Path to the Seurat object (.Rds file) and 2. Clustering method either seurat or monocle")
}
print(args)
brain <- readRDS(args[1])
clusterMethod <- args[2]

clusterColumn <- paste0(clusterMethod,"_clusters")
sampleName    <- gsub("Seurat_","",basename(tools::file_path_sans_ext(args[1])))


setwd("clustering_snakemake/01.clustering/")
dir.create(paste0("bigwig/",sampleName),recursive=TRUE)
setwd(paste0("bigwig/",sampleName))


if (!clusterMethod %in% c("seurat","monocle")){
  stop("Cluster method must be either seurat or monocle ")
}

if (grepl("_UMI",args[1])){
  fragments <- rtracklayer::import("../../../../data/combined/outs/fragments.tsv.gz",format = "bed")
} else if (grepl("_noUMI",args[1])) {
  fragments <- rtracklayer::import("../../../../data/combined/outs/noUMI/fragments.tsv.gz",format = "bed")
} else {
  stop("Can't determine whether UMI condition (_UMI or _noUMI in .Rds filename) ")
}

# Merge metadata
barcode.export <- brain@meta.data
rownames(barcode.export) <- paste(barcode.export$sample,barcode.export$barcode,sep="_")


# Filter metadata to clustered cells
barcode.export <- barcode.export[!is.na(barcode.export[,clusterColumn]),]
barcode.export[,clusterColumn] <- as.character(barcode.export[,clusterColumn])

# Split metadata by cluster
barcodes.by.cluster <-lapply(sort(unique(barcode.export[,clusterColumn])),function(x){
  rownames(barcode.export[barcode.export[,clusterColumn] == x,])
})

# Get reads per cluster from fragments
fragments.by.cluster <- lapply(barcodes.by.cluster,function(x){
  fragments[fragments$name %in% x]
})

# Convert reads to coverage
coverage.by.cluster <- lapply(fragments.by.cluster,GenomicRanges::coverage)

# Normalize
coverage.by.cluster <- lapply(seq(coverage.by.cluster),function(x){
  coverage.by.cluster[[x]] <- coverage.by.cluster[[x]]/(length(fragments.by.cluster[[x]]) / 10^6)
  coverage.by.cluster[[x]]
})

# Export
lapply(seq(coverage.by.cluster),function(x){
  rtracklayer::export.bw(object = coverage.by.cluster[[x]],con = paste0(clusterMethod,"_cluster_",x,".bw"))
})

fragments.all.clusters <- do.call("c",fragments.by.cluster)
coverage.all.clusters  <- coverage(fragments.all.clusters)
coverage.all.clusters  <- coverage.all.clusters/length(fragments.all.clusters)

rtracklayer::export(object = coverage.all.clusters,con = paste0(clusterMethod,"_all_clusters.bw"),format = "bw")

############### Export top 200 cells per cluster

# dir.create("individual_cells")
# setwd("individual_cells")
# 
# 
# clusters_top200 <- lapply(sort(unique(barcode.export[,clusterColumn])),function(x,n=100){
#   cluster.meta <- barcode.export[barcode.export[,clusterColumn] == x,]
#   cluster.meta <- cluster.meta[order(cluster.meta$all_unique_MB,decreasing=TRUE),]
#   cluster.meta <- head(cluster.meta,n)
#   cluster.meta
# })
# 
# lapply(seq(clusters_top200),function(x){
#   fragments.cluster <- fragments[fragments$name %in% rownames(clusters_top200[[x]])]
#   fragments.cluster.ls <- split(fragments.cluster,fragments.cluster$name)
#   fragments.cluster.ls.coverage <- lapply(fragments.cluster.ls,coverage)
#   lapply(seq(fragments.cluster.ls.coverage),function(y){
#     rtracklayer::export(object=fragments.cluster.ls.coverage[[y]],con=paste0("Cluster_",x,"_",names(fragments.cluster.ls.coverage)[y],"_reads.bw"))
#     y
#     })
# })



