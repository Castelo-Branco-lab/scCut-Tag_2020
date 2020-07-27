library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(viridis)
library(ggplot2)
library(monocle3)


# Debug
#args <- c("clustering_snakemake/01.clustering/Seurat_5000_UMI.Rds","monocle")
args <- commandArgs(trailingOnly = TRUE)

print(args)
if (length(args) != 2){
  stop("Must provide 2 arguments: 1. Path to the Seurat object (.Rds file) and 2. Clustering method either seurat or monocle")
}

brain <- readRDS(args[1])
clusterMethod <- args[2]

clusterColumn <- paste0(clusterMethod,"_clusters")
sampleName    <- gsub("Seurat_","",basename(tools::file_path_sans_ext(args[1])))


setwd("clustering_snakemake/01.clustering/")
dir.create(paste0("markers/",sampleName),recursive=TRUE)
setwd(paste0("markers/",sampleName))


######################################################################
#################### Find markers for populations ####################
######################################################################
brain@active.ident <- brain@meta.data[,paste0(clusterMethod,"_clusters")]
names(brain@active.ident) <- names(brain$barcode)

# Filter away clusters with < 200 cells
clusters_200cells <- names(table(brain@active.ident)[table(brain@active.ident) > 200])
brain <- brain[,brain@active.ident %in% clusters_200cells]

gene.coords <- ensembldb::genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
lncRNA.coords <- ensembldb::genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "lincRNA")
gene.coords <- c(gene.coords,lncRNA.coords)

seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

# Flatten the overlapping genes and extend by 2kb upstream of promoters
genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat,upstream = 2000)

# Retrieve gene names from the original annotation (lost because of flatenning)
genebodyandpromoter.coords.flat$name<- gene.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$gene_name


####################
brain <- ScaleData(object = brain, features = rownames(brain),assay = "peaksMB")
diffPeaks <- FindAllMarkers(brain,assay = "peaksMB",min.pct = 0.025,logfc.threshold = 0.2)
diffPeaks$closest_gene <- ClosestFeature(diffPeaks$gene,annotation = genebodyandpromoter.coords.flat)$name

topMarkers <- diffPeaks %>% dplyr::filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv2(x = topMarkers, paste0(clusterMethod,"_top_markers_peaks.tsv"),quote = FALSE,row.names = FALSE)

topMarkers <- diffPeaks %>% dplyr::filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

# Remove duplicate markers - keeps one instance of the marker only
topMarkers <- topMarkers[!duplicated(topMarkers$gene),]


#brain@active.ident <- factor(brain@active.ident,levels=c(0,2,1,3,4,5))
p1 <- DoHeatmap(obj = brain, features = unique(topMarkers$gene), assay = "peaksMB",raster = FALSE) + 
  scale_fill_viridis_c() 

heatmap.ylabs <- ggplot_build(p1)$layout$panel_params[[1]]$y.labels
topMarkers <- as.data.frame(topMarkers)
rownames(topMarkers) <- topMarkers$gene
topMarkers <- topMarkers[heatmap.ylabs,]

p1 + scale_y_discrete(labels=paste(topMarkers$gene,topMarkers$closest_gene))

ggsave(paste0(clusterMethod,"_heatmap_peaks_top.pdf"),width=12,height=12)
ggsave(paste0(clusterMethod,"_heatmap_peaks_top.png"),width=12,height=12)


########################################
######### Gene activity matrix #########
########################################

if (grepl("_UMI",args[1])){
  fragments_file <- "../../../../data/combined/outs/fragments.tsv.gz"
} else if (grepl("_noUMI",args[1])) {
  fragments_file <- "../../../../data/combined/outs/noUMI/fragments.tsv.gz"
} else {
  stop("Can't determine whether UMI condition (_UMI or _noUMI in .Rds filename) ")
}

gene.activities <- FeatureMatrix(
  fragments = fragments_file,
  features = genebodyandpromoter.coords.flat,
  cells = colnames(brain),
  chunk = 10
)


gene.key <- genebodyandpromoter.coords.flat$name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords.flat)
rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
gene.activities <- gene.activities[rownames(gene.activities)!="",]

brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)

###########################################
###### Find markers and plot heatmap ######
###########################################

DefaultAssay(brain) <- 'RNA'
brain <- ScaleData(object = brain, features = rownames(brain),assay = "RNA",)

diffGenes <- FindAllMarkers(brain,assay = "RNA",min.pct = 0.05,logfc.threshold = 0.1)
topMarkers <- diffPeaks %>% dplyr::filter(avg_logFC > 0) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv2(x = topMarkers, paste0(clusterMethod,"_top_markers_genes.tsv"),quote = FALSE,row.names = FALSE)



topMarkers <- diffGenes %>% group_by(cluster) %>% dplyr::filter(avg_logFC > 0) %>% top_n(n = 25, wt = avg_logFC)
topMarkers <- topMarkers[!duplicated(topMarkers$gene),]


DoHeatmap(obj = brain,topMarkers$gene) + scale_fill_viridis_c()
ggsave(paste0(clusterMethod,"_heatmap_genes_top.pdf"),width=10,height=10)
ggsave(paste0(clusterMethod,"_heatmap_genes_top.png"),width=10,height=10)


###################################
###### Monocle3 gene modules ######
###################################
if(clusterMethod == "monocle"){

  monocle_path <- paste0("../../",basename(gsub("Seurat","Monocle",args[1])))
  cds <- readRDS(monocle_path)
  
  gene_module_df <- find_gene_modules(cds, resolution=0.01,max_components=8)
  gene_module_df$gene_name <- ClosestFeature(gene_module_df$id,annotation=genebodyandpromoter.coords.flat)$name
  
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=partitions(cds))
  
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
  
  filter_modules_heatmap <- function(agg_mat,nmodules) {
    filter_modules = names(head(sort(apply(agg_mat,1,function(x){max(x)-min(x) }),decreasing=TRUE),nmodules))
    agg_mat_filtered <- agg_mat[filter_modules,]
    
    pheatmap::pheatmap(agg_mat_filtered, cluster_rows=TRUE, cluster_cols=TRUE,
                       scale="row",  clustering_method="ward.D2",
                       fontsize=6,filename=paste0("Heatmap_modules_",nmodules,".pdf"))
  }
    
  filter_modules_heatmap(agg_mat,10)
  filter_modules_heatmap(agg_mat,20)
  filter_modules_heatmap(agg_mat,30)
  filter_modules_heatmap(agg_mat,40)
  filter_modules_heatmap(agg_mat,50)
  filter_modules_heatmap(agg_mat,100)
  filter_modules_heatmap(agg_mat,200)
  filter_modules_heatmap(agg_mat,500)
  
  genes_per_module <- aggregate(gene_module_df$gene_name,by=list(gene_module_df$module),FUN="paste")
  genes_per_module$x <- lapply(genes_per_module$x,unique)
  genes_per_module$x <- unlist(lapply(genes_per_module$x,paste,collapse=", "))
  colnames(genes_per_module) <- c("module number","genes")
  
  write.csv2(x = gene_module_df,file = "Modules_annotations_long.csv",quote = FALSE,row.names = FALSE)
  write.csv2(x=genes_per_module,file="Modules_annotations_aggregated.csv",quote = FALSE,row.names = FALSE)

}









