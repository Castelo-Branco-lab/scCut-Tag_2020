library(argparse)
library(Seurat)
library(Signac)
#library(EnsDb.Mmusculus.v79);
#library(viridis);
#library(ggplot2);
#library(dplyr);
#library(BSgenome.Mmusculus.UCSC.mm10)

###############################
parser <- ArgumentParser()

parser$add_argument("-s", "--sample", type="character", default='foo', 
                    help="sample name [as in config file key]")

parser$add_argument("-m", "--metadata", type="character", default='foo', 
                    help="path to the pre-processed metadata [output from pick_cells.R]")

parser$add_argument("-c", "--config", type="character", default='config/config.yaml', 
                    help="maximum number of reads in cell")

parser$add_argument("-o", "--out_prefix", type="character", default="10000", 
                    help="folder for the output in clustering_snakemake folder")

parser$add_argument("-w", "--window", type="numeric", default="10000", 
                    help="width of a window")

args <- parser$parse_args()

############################ Filter the dataset

metadata <- read.csv(file = paste0(args$out_prefix,'metadata.csv'))
metadata <- metadata[metadata$passedMB,]
rownames(metadata) <- metadata$barcode

############################ Create peak counts matrix
cat("*** Creating bins/peaks matrices \n")

peaks_file = paste0('results/',args$sample,'/macs/broad/',args$sample,'_peaks.broadPeak')
if (!file.exists(peaks_file)) {stop(paste0("Peaks file does not exist:: ", peaks_file))}
peaks <- rtracklayer::import(peaks_file)

counts.matrix.peaks <- FeatureMatrix(fragments = fragments,
                                     features = peaks,
                                     cells = metadata$barcode,
                                     chunk = 50)


counts.matrix.bins <- GenomeBinMatrix(fragments = fragments,
                                      genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
                                      binsize = window,
                                      cells = metadata$barcode)

########################## Create Seurat object
min_features = 1
min_cells    = 1

seurat_object <- CreateSeuratObject(counts = counts.matrix.bins,
                     project = args$sample,
                     assay = paste0('bins_',args$window),
                     meta.data = metadata,
                     min.features = min_features,
                     min.cells = min_cells)


seurat_object[['peaks']] <- CreateAssayObject(counts = counts.matrix.peaks[,colnames(counts.matrix.peaks) %in% colnames(seurat_object)])

# Filter blacklist cells
seurat_object$blacklist_ratio <- seurat_object$blacklist_region_fragments / seurat_object$all_unique_MB
seurat_object                 <- seurat_object[,seurat_object$blacklist_region_fragments < 5]

# Add sample id to cell names
seurat_object            <- RenameCells(object = seurat_object,add.cell.id = args$sample)

seurat_object$antibody   <- config$sample[[args$sample]]$Antibody
seurat_object$GFP        <- config$sample[[args$sample]]$GFP
seurat_object$Age        <- config$sample[[args$sample]]$Age


########## Create Gene activity matrix
cat("*** Create gene matrix for mm10 \n")
gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat,
                                 cells = gsub(paste0(args$sample,"_"),"",colnames(seurat_object)))

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)
rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]

gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]
colnames(gene.matrix) <- paste0(args$sample,"_",colnames(gene.matrix))

seurat_object[['GA']] <- CreateAssayObject(counts = gene.matrix)

######### Save the object
cat("*** Save the object \n")
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

########################## Try Lustering using Seurat
cat("*** Clustering and dimensionality reduction \n")

DefaultAssay(seurat_object) <- paste0('bins_',args$window)

# The dimensionality reduction might fail, especially for low complexity datasets, so put the code into try to still save the objects
try({
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object,min.cutoff = 'q20')
seurat_object <- RunSVD(
  seurat_object,
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
  )


seurat_object <- RunUMAP(seurat_object, dims = 2:ndim, reduction = 'lsi')

seurat_object <- FindNeighbors(
  object = seurat_object,
  reduction = 'lsi',
  dims = 2:ndim
)


seurat_object <- FindClusters(
  object = seurat_object,
  algorithm = "leiden",
  resolution = 0.3,
  verbose = TRUE
)


p1 <- DimPlot(seurat_object,group.by = 'ident')
p2 <- FeaturePlot(seurat_object,'blacklist_region_fragments')
p3 <- FeaturePlot(seurat_object,'logUMI')



ggsave(plot= p1 + p2 + p3,
       filename = paste0(args$out_prefix,'Seurat_clustering.png'),
       width=15,height=5)
})

######### Save the final object
cat("*** Save the object \n")
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

