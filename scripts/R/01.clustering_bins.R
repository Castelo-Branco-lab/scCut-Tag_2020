# Final clustering for H3K4me3 dataset
  
  library(argparse)
  library(monocle3)
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(viridis)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(BSgenome.Mmusculus.UCSC.mm10)
  
  set.seed(1234)
  rm(list=ls())
  
  
  ########### Arguments parser
    
    parser <- ArgumentParser()
    
    parser$add_argument("-l", "--reads_min", type="double", default=3.5, 
                        help="minimum number of reads in cell [log10 scale]")
    
    parser$add_argument("-m", "--reads_max", type="double", default=5.5, 
                        help="maximum number of reads in cell [log10 scale]")
    
    parser$add_argument("-q", "--peaks_min", type="double", default=0.215, 
                        help="maximum number of reads in cell")
    
    parser$add_argument("-r", "--peaks_max", type="double", default=0.95, 
                        help="maximum number of reads in cell")
    
    parser$add_argument("-d", "--ndim", type="integer", default=8, 
                        help="number of dimensions for UMAP reduction")
    
    parser$add_argument("-c", "--res_cluster", type="double", default=0.2, 
                        help="clustering resolution")
    
    parser$add_argument("-w", "--bin_width", type="integer", default=10000, 
                        help="width of the bin")
    
    
    
    
    args <- parser$parse_args()
    print(args)
    
    cutoff_reads_min            = 10^args$reads_min
    cutoff_reads_max            = 10^args$reads_max
    cutoff_peak_percentage_low  = args$peaks_min
    cutoff_peak_percentage_high = args$peaks_max
    
    ndim = args$ndim
    res_clustering = args$res_cluster
    window                      = args$bin_width
    
    assay = "peaksMB"
    
  
##############################################################
dir.create("clustering_snakemake")
setwd("clustering_snakemake")

dir.create("01.clustering_bins")
setwd("01.clustering_bins")

dir.create("figures")
###############################################################

if(grepl("monod.*.mbb.ki.se", Sys.info()["nodename"])) {
  Sys.setenv(PATH=paste0('/home/marek/anaconda3/envs/CT/bin:', Sys.getenv('PATH')))
  Sys.setenv(RETICULATE_PYTHON = "/home/marek/anaconda3/envs/CT/bin/python")
}

if(Sys.info()["nodename"] == "KI-DGKT5HVTGG7F"){
  Sys.setenv(PATH=paste0('/usr/local/anaconda3/bin:', Sys.getenv('PATH')))
  Sys.setenv(RETICULATE_PYTHON = "/usr/local/anaconda3/bin/python3")
}
###############################################################

#### Read config file
config <- yaml::read_yaml("../../config.yaml")
samples <- names(config$samples)
files   <- unlist(config$samples)
antibody <- strsplit(samples,"_N")[[1]][1]

if(!grepl("monod.*.mbb.ki.se",Sys.info()['nodename'])){
  files <- gsub("/data/proj/GCB_MB","~/mount",files)
}

########## Filter the barcodes

all_barcodes  <- paste0("../../barcode_statistics/",samples,"/all_barcodes_broad.txt")
peak_barcodes <- paste0("../../barcode_statistics/",samples,"/peaks_barcodes_broad.txt")
md_files      <- paste0(files,"/outs/singlecell.csv")


barcode.ls.10X = lapply(seq(all_barcodes), function(i){
  barcodes = read.csv(
    md_files[i], 
    head=TRUE
  )
  barcodes = barcodes[2:nrow(barcodes),]
  #barcodes = barcodes[barcodes$cell_id != "None",]          # Filter away cell barcodes that were not called by cellranger
  barcodes$logUMI = log10(barcodes$passed_filters + 1)
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1)
  barcodes$peak_region_ratio = (barcodes$peak_region_fragments+1) / (barcodes$passed_filters + 1)
  barcodes
})

barcode.ls.10X <- lapply(1:length(barcode.ls.10X),function(x) {
  all_barcodes <- read.table(file=all_barcodes[x])
  peak_barcodes <- read.table(file=peak_barcodes[x])
  bcd <- merge(all_barcodes,peak_barcodes,by="V2")  
  colnames(bcd) <- c("barcode","all_unique_MB","peak_MB")
  bcd$peak_ratio_MB <- bcd$peak_MB/bcd$all_unique_MB
  bcd$sample <- samples[x]
  return(merge(barcode.ls.10X[[x]],bcd,by.x="barcode",by.y="barcode"))
})

barcode.ls.10X <- lapply(barcode.ls.10X,function(x){
  x$is__cell_barcode <- as.factor(x$is__cell_barcode)
  x
})
################ MB filtering

barcode.ls.10X <- lapply(barcode.ls.10X,function(x){
  x[,"passedMB"] <- FALSE
  x[x$all_unique_MB > cutoff_reads_min &
      x$all_unique_MB < cutoff_reads_max &
      x$peak_ratio_MB > cutoff_peak_percentage_low &
      x$peak_ratio_MB < cutoff_peak_percentage_high,"passedMB"] <- TRUE
  x
})

################ Cell picking scatterplot nreads ~ percent in peaks


p1 <- lapply(seq(barcode.ls.10X),function(x){
  ggplot(data=barcode.ls.10X[[x]]) + 
    geom_point(aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=is__cell_barcode),size=0.1) + 
    scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(barcode.ls.10X[[x]]$passedMB)),"FALSE")) +
    theme(legend.position="bottom")    
})

p2 <- lapply(seq(barcode.ls.10X),function(x){
  ggplot(data=barcode.ls.10X[[x]]) + 
    geom_point(aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=is__cell_barcode),size=0.1) +
    scale_color_manual(values=c("black","gold")) +
    theme(legend.position="bottom")    
  
})


plots <- do.call('c',list(p1,p2))
g <- do.call("grid.arrange", c(plots,ncol=floor(sqrt(length(plots)))))
ggsave(filename = paste0("figures/Bins_",window,"_scatterplot_cells_picking_10x.pdf"),plot = g,width=10,heigh=10)

#

p1 <- lapply(seq(barcode.ls.10X),function(x){
  ggplot(data=barcode.ls.10X[[x]]) + 
    geom_point(aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=passedMB),size=0.1) + 
    geom_hline(yintercept = c(cutoff_peak_percentage_high,cutoff_peak_percentage_low)) + 
    geom_vline(xintercept = log10(c(cutoff_reads_min,cutoff_reads_max))) + 
    scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(barcode.ls.10X[[x]]$passedMB)),
                                                         "FALSE")) +
    theme(legend.position="bottom")    
})

p2 <- lapply(seq(barcode.ls.10X),function(x){
  ggplot(data=barcode.ls.10X[[x]]) + 
    geom_point(aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=passedMB),size=0.1) +
    scale_color_manual(values=c("black","gold")) +
    theme(legend.position="bottom")    
  
})


plots <- do.call('c',list(p1,p2))
g <- do.call("grid.arrange", c(plots,ncol=floor(sqrt(length(plots)))))
ggsave(filename = paste0("figures/Bins_",window,"_scatterplot_cells_picking_custom.pdf"),plot = g,width=10,heigh=10)


rm(list = c("g","p1","p2","plots"))
############################ Filter the dataset

barcode.ls.10X.filter <- lapply(barcode.ls.10X,function(x){
  x <- x[x$all_unique_MB > cutoff_reads_min &
           x$all_unique_MB < cutoff_reads_max &
           x$peak_ratio_MB > cutoff_peak_percentage_low &
           x$peak_ratio_MB < cutoff_peak_percentage_high,]
  rownames(x) <- x$barcode
  x
})

############################ Merge barcode metadata
barcode.10X.filter <- do.call('rbind',barcode.ls.10X.filter)
rownames(barcode.10X.filter) <- paste0(barcode.10X.filter$sample,"_",barcode.10X.filter$barcode)


############################ Create peak activities matrix

peaks <- rtracklayer::import("../../macs2_peaks.combined/broad/combined_peaks.broadPeak")
fragments <- "../../data/combined/outs/fragments.tsv.gz"


bin.activities <- GenomeBinMatrix(fragments = fragments, 
                             genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
                             binsize = window)


# Filter bins, remove bottom 5% and top 2%
#min_features = quantile(Matrix::rowSums(bin.activities),0.05)
#max_features <- quantile(Matrix::rowSums(bin.activities),0.98)


bin.activities.filtered <- bin.activities[,which(colnames(bin.activities) %in% rownames(barcode.10X.filter))] # Filter cells that pass in standard Seurat clustering

#bin.activities.filtered <- bin.activities.filtered[Matrix::rowSums(bin.activities.filtered) > min_features,] # min features filtering 
#bin.activities.filtered <- bin.activities.filtered[Matrix::rowSums(bin.activities.filtered) < max_features,] # min features filtering 
#dim(bin.activities.filtered)

########################## Create Seurat object 
min_features = 1
min_cells    = 1

brain <- CreateSeuratObject(counts = bin.activities.filtered,
                            project = antibody,
                            assay = assay,
                            meta.data = barcode.10X.filter,
                            min.features = min_features,
                            min.cells = min_cells)



########################## CLuster using Seurat


brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain,min.cutoff = 'q5')
brain <- RunSVD(
  brain,
  reduction.key = 'LSI_',
  reduction.name = 'lsi', 
  n = 50
)


brain <- RunUMAP(brain, dims = 2:ndim, reduction = 'lsi')

brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:ndim
)


brain <- FindClusters(
  object = brain,
  algorithm = "leiden",
  resolution = res_clustering,
  verbose = TRUE
)

p1 <- DimPlot(brain,group.by = 'sample')
p2 <- DimPlot(brain,group.by = 'ident')

ggsave(plot=p1+p2,filename = paste0("figures/01.Bins_",window,"_Seurat_clustering.pdf"),width=12,height=6)



########################## Cluster using Monocle3

ndim = 8
plot.group = "sample"

cds <- new_cell_data_set(expression_data = brain[[assay]]@counts, cell_metadata = brain@meta.data)
cds <- preprocess_cds(cds, num_dim = ndim, method = "LSI")
cds <- reduce_dimension(cds,preprocess_method = "LSI", reduction_method="UMAP")
cds <- cluster_cells(cds)
cds <- invisible(learn_graph(cds))

p1 <- plot_cells(cds)
p2 <- plot_cells(cds,color_cells_by = "sample")

ggsave(plot=p1+p2, filename=paste0("figures/01.Bins_",window,"_Monocle_clustering.pdf"),device='pdf',width=12,height=6)

brain[['monocle_lsi']] <- CreateDimReducObject(embeddings=reducedDims(cds)[['LSI']],assay='peaksMB',key='monoclelsi_')
brain[['monocle_umap']] <- CreateDimReducObject(embeddings=reducedDims(cds)[['UMAP']],assay='peaksMB',key='monocleumap_')


monocle.clusters <- cds@clusters[["UMAP"]]$clusters
brain$monocle_clusters <- monocle.clusters[colnames(brain)]


# 
p1 <- DimPlot(brain,group.by = 'seurat_clusters')
p2 <- DimPlot(brain,group.by = 'monocle_clusters')

ggsave(plot=p1+p2, filename = paste0("figures/01.Bins_",window,"_Seurat_monocle_in_monocle_dims.pdf"),width=12,height=6,device='pdf')

#####################
source("~/bin/snaptools_scripts/snakemake/R/01.Clustering_functions.R")

plots <- lapply(c(3:15,20,25,30),function(x){
  clusterSeurat(brain,x)
})

g <- do.call("grid.arrange",plots)
ggsave(filename = paste0("01.Bins_Seurat_",window,"_ndim_matrix.pdf"),plot = g,width=12,heigh=12)


plots <- lapply(c(3:15,20,25,30),function(x){
  clusterMonocle(brain,"peaksMB",x)
}) 

g <- do.call("grid.arrange",plots)
ggsave(filename = paste0("01.Bins_Monocle_",window,"_ndim_matrix.pdf"),plot = g,width=12,heigh=12)
#####################


saveRDS(brain,paste0("01.Seurat_object_bins_",window,".Rds"))
