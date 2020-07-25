# Final clustering for H3K4me3 dataset
suppressMessages({
library(argparse);
library(Seurat);
library(Signac);
library(EnsDb.Mmusculus.v79);
library(EnsDb.Hsapiens.v86);
library(viridis);
library(ggplot2);
library(dplyr);
library(gridExtra);
library(BSgenome.Mmusculus.UCSC.mm10);}
)

set.seed(1234)
rm(list=ls())


########### Arguments parser

parser <- ArgumentParser()

parser$add_argument("-l", "--reads_min", type="double", default=4, 
                    help="minimum number of reads in cell [log10 scale]")

parser$add_argument("-m", "--reads_max", type="double", default=5.5, 
                    help="maximum number of reads in cell [log10 scale]")

parser$add_argument("-q", "--peaks_min", type="double", default=0.25, 
                    help="maximum number of reads in cell")

parser$add_argument("-r", "--peaks_max", type="double", default=0.7, 
                    help="maximum number of reads in cell")

parser$add_argument("-d", "--ndim", type="integer", default=30, 
                    help="maximum number of reads in cell")

parser$add_argument("-c", "--res_cluster", type="double", default=0.2, 
                    help="maximum number of reads in cell")

parser$add_argument("-o", "--out_prefix", type="character", default="100000_UMI", 
                    help="folder for the output in clustering_snakemake folder")

parser$add_argument("-t", "--feature", type="character", default="bins", 
                    help="type of features choose one of: [bins,peaks]")

parser$add_argument("-g", "--genome", type="character", default="mm10", 
                    help="genome either mm10 or GRCh38")

parser$add_argument("-p", "--peaks_file", type="character", default="macs2_peaks.combined/broad/combined_peaks.broadPeak", 
                    help="path to the peaks file [.bed]")

parser$add_argument("-f", "--fragments_file", type="character", default="data/combined/outs/fragments.tsv.gz", 
                    help="path to the fragments file")

parser$add_argument("-w", "--window", type="integer", default=100000, 
                    help="width of the genome window (if -t bins)")

parser$add_argument("-1", "--cells_selection", type="logical", default=FALSE, 
                    help="Quit after plotting the cell depth-percentage in peaks plot to adjust the parameters")

args <- parser$parse_args()
print(args)

cutoff_reads_min            = 10^args$reads_min
cutoff_reads_max            = 10^args$reads_max
cutoff_peak_percentage_low  = args$peaks_min
cutoff_peak_percentage_high = args$peaks_max

ndim = args$ndim
res_clustering = args$res_cluster
window = args$window

fragments <- args$fragments_file
assay = "peaksMB"


##############################################################
dir.create("clustering_snakemake/01.clustering/figures",recursive=TRUE)

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

#### Read annotation

if (args$genome == 'mm10')  {ensdb = EnsDb.Mmusculus.v79}
if (args$genome == 'GRCh38'){ensdb = EnsDb.Hsapiens.v86}

gene.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "protein_coding")
lncRNA.coords <- ensembldb::genes(ensdb, filter = ~ gene_biotype == "lincRNA")
gene.coords <- c(gene.coords,lncRNA.coords)

seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

# Flatten the overlapping genes and extend by 2kb upstream of promoters
genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat,upstream = 2000)

# Retrieve gene names from the original annotation (lost because of flatenning)
genebodyandpromoter.coords.flat$name<- gene.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$gene_name



#### Read config file
if (!file.exists("config.yaml")) {
  stop("You are in the wrong the directory !")
}

config <- yaml::read_yaml("config.yaml")
samples <- names(config$samples)
files   <- unlist(config$samples)
antibody <- strsplit(samples,"_N")[[1]][1]

if(!grepl("monod.*.mbb.ki.se",Sys.info()['nodename'])){
  files <- gsub("/data/proj/GCB_MB","~/mount",files)
}

########## Filter the barcodes

all_barcodes  <- paste0("barcode_statistics/",samples,"/all_barcodes_broad.txt")
peak_barcodes <- paste0("barcode_statistics/",samples,"/peaks_barcodes_broad.txt")
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

if(!file.exists("clustering_snakemake/01.clustering/figures/scatterplot_cells_picking_10x.pdf")){
  p1 <- lapply(seq(barcode.ls.10X),function(x){
    ggplot(data=barcode.ls.10X[[x]]) + 
      geom_point(aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=is__cell_barcode),size=0.1) + 
      scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(barcode.ls.10X[[x]]$passedMB)),
                                                             "FALSE")) +
      theme(legend.position="bottom")    
  })
  
  p2 <- lapply(seq(barcode.ls.10X),function(x){
    ggplot(data=barcode.ls.10X[[x]]) + 
      geom_point(aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=is__cell_barcode),size=0.1) +
      scale_color_manual(values=c("black","gold"))
  })
  
  
  plots <- do.call('c',list(p1,p2))
  g <- do.call("grid.arrange", c(plots,ncol=floor(sqrt(length(plots)))))
  ggsave(filename = "clustering_snakemake/01.clustering/figures/scatterplot_cells_picking_10x.pdf",
         plot = g,
         width=10,heigh=10)

}
#
if(!file.exists("clustering_snakemake/01.clustering/figures/scatterplot_cells_picking_custom.pdf")){
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
      scale_color_manual(values=c("black","gold"))
  })
  
  
  plots <- do.call('c',list(p1,p2))
  g <- do.call("grid.arrange", c(plots,ncol=floor(sqrt(length(plots)))))
  ggsave(filename = "clustering_snakemake/01.clustering/figures/scatterplot_cells_picking_custom.pdf",
         plot = g,
         width=10,heigh=10)
  rm(list = c("g","p1","p2","plots"))
}

################# Export bw selected / unselected
fragments_gr <- rtracklayer::import(fragments,format = "bed")

lapply(barcode.ls.10X,function(x){
  x.barcode_pass <- paste0(x$sample,"_",x$barcode)[x$passedMB]
  x.barcode_nopass <- paste0(x$sample,"_",x$barcode)[!x$passedMB]
  
  print("Retrieving reads")
  fragments.pass <- fragments_gr[fragments_gr$name %in% x.barcode_pass]
  fragments.nopass <- fragments_gr[fragments_gr$name %in% x.barcode_nopass]
  
  print(c(length(fragments.pass),length(fragments.nopass)))
  
  print("Calculating coverage")
  coverage.pass <- GenomicRanges::coverage(fragments.pass)
  coverage.nopass <- GenomicRanges::coverage(fragments.nopass)
  
  print("Normalizing")
  coverage.pass <- coverage.pass/length(fragments.pass)
  coverage.nopass <- coverage.nopass/length(fragments.nopass)
  
  print("Exporting")
  #rtracklayer::export(object=coverage.pass,con = paste0("~/snakemake/H3K4me3/01.clustering/bigwig/",x$sample[1],"_pass.bw"))
  #rtracklayer::export(object=coverage.nopass,con = paste0("~/snakemake/H3K4me3/01.clustering/bigwig/",x$sample[1],"_nopass.bw"))
  
  
  rtracklayer::export(object=coverage.pass,con = paste0("bigwig/",x$sample[1],"filter_pass.bw"))
  rtracklayer::export(object=coverage.nopass,con = paste0("bigwig/",x$sample[1],"filter_nopass.bw"))
  
  
})


############ Stop here if cell_selection is selected
if(args$cells_selection){
  quit("no",0)
}

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


############################ Create peak counts matrix

if (!file.exists(fragments)) {stop(paste0("Fragments file does not exist: ",fragments))}

if(args$feature == "peaks"){
  if (!file.exists(args$peaks_file)) {stop(paste0("Fragments file does not exist:: ",args$peaks_file))}
  peaks <- rtracklayer::import(args$peaks_file)
  
  counts.matrix <- FeatureMatrix(
    fragments = fragments,
    features = peaks,
    cells = rownames(barcode.10X.filter),
    chunk = 50
  )
}



if(args$feature == "bins"){
    counts.matrix <- GenomeBinMatrix(fragments = fragments, 
                                      genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
                                      binsize = window)
  }



counts.matrix <- counts.matrix[,colnames(counts.matrix) %in% rownames(barcode.10X.filter)]

########################## Create Seurat object 
min_features = 1
min_cells    = 1

brain <- CreateSeuratObject(counts = counts.matrix,
                     project = antibody,
                     assay = assay,
                     meta.data = barcode.10X.filter,
                     min.features = min_features,
                     min.cells = min_cells)

brain <- brain[,Matrix::colSums(brain[[assay]]@counts) >=min_cells]
brain <- brain[Matrix::rowSums(brain[[assay]]@counts,) >=min_features]


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

ggsave(plot=p1+p2,
       filename = paste0("clustering_snakemake/01.clustering/figures/","Seurat_ndim_",ndim,"_",args$out_prefix,".pdf"),
       width=12,height=6)

########## Create Gene activity matrix 

gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat,cells = colnames(brain))

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)
rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]

gene.matrix   <- gene.matrix[rownames(gene.matrix) != "",]
brain[['GA']] <- CreateAssayObject(counts = gene.matrix)


########################## Cluster using Monocle3
# gene.metadata <- data.frame(gene_short_name = rownames(brain))
# rownames(gene.metadata) <- rownames(brain)
#
# cds <- new_cell_data_set(expression_data = brain[[assay]]@counts, cell_metadata = brain@meta.data,gene_metadata = gene.metadata)
# cds <- preprocess_cds(cds, num_dim = ndim, method = "LSI")
# cds <- reduce_dimension(cds,preprocess_method = "LSI", reduction_method="UMAP")
# cds <- cluster_cells(cds)
# cds <- invisible(learn_graph(cds))
#
# p1 <- plot_cells(cds)
# p2 <- plot_cells(cds,color_cells_by = "sample")
#
# ggsave(plot=p1+p2,
#        filename = paste0("clustering_snakemake/01.clustering/figures/","Monocle_ndim_",ndim,"_",args$out_prefix,".pdf"),
#        width=12,height=6)
#
# brain[['monocle_lsi']] <- CreateDimReducObject(embeddings=reducedDims(cds)[['LSI']],assay='peaksMB',key='monoclelsi_')
# brain[['monocle_umap']] <- CreateDimReducObject(embeddings=reducedDims(cds)[['UMAP']],assay='peaksMB',key='monocleumap_')
#
#
# monocle.clusters <- cds@clusters[["UMAP"]]$clusters
# brain$monocle_clusters <- monocle.clusters[colnames(brain)]
#
#
# #
# p1 <- DimPlot(brain,group.by = 'seurat_clusters')
# p2 <- DimPlot(brain,group.by = 'monocle_clusters')
#
# ggsave(plot=p1+p2, filename=paste0("clustering_snakemake/01.clustering/figures/Monocle_ndim_",ndim,"_",args$out_prefix,".pdf"),width=12,height=6,device='pdf')

# #####################
# source("~/bin/snaptools_scripts/snakemake/R/01.Clustering_functions.R")
#
# plots <- lapply(c(3:15,20,25,30),function(x){
#   clusterSeurat(brain,x)
#   })
#
# g <- do.call("grid.arrange",plots)
# ggsave(filename = paste0("clustering_snakemake/01.clustering/figures/","Seurat_matrix_",args$out_prefix,".pdf"),
#        plot = g,
#        width=12,heigh=12)
#
#
# plots <- lapply(c(3:15,20,25,30),function(x){
#   clusterMonocle(brain,"peaksMB",x)
# })
#
# g <- do.call("grid.arrange",plots)
# ggsave(filename = paste0("clustering_snakemake/01.clustering/figures/","Monocle_matrix_",args$out_prefix,".pdf"),
#        plot = g,
#        width=12,heigh=12)

####################
# saveRDS(cds,
#         paste0("clustering_snakemake/01.clustering/Monocle_",args$out_prefix,".Rds")
# )

saveRDS(brain,
        paste0("clustering_snakemake/01.clustering/Seurat_",args$out_prefix,".Rds")
)

