library(argparse)
library(Seurat)
library(Signac)
library(GenomicFeatures)
library(yaml)
library(rtracklayer)
library(viridis);
library(ggplot2);
library(funr)


source(paste0(dirname(funr::sys.script()),"/func.R"))
ndim = 40

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

parser$add_argument("-w", "--window", type="integer", default="10000", 
                    help="width of a window")

args <- parser$parse_args()

####### Load the config
config    <- yaml::yaml.load_file(input = args$config)

# Get genome version
genome_version <- config$samples[[args$sample]]$Version
if (is.null(genome_version)){
  genome_version <- config$general$Version
}

if (is.null(genome_version)) {
  cat("*** Error: genome version needs to be specified in config.yaml either sample specific or in general ***")
}

############################ Filter the dataset

metadata <- read.csv(file = args$metadata,stringsAsFactors = FALSE)
metadata <- metadata[metadata$passedMB,]
rownames(metadata) <- metadata$barcode

######## Fragments
fragments.path <- paste0(config$samples[[args$sample]]$cellranger_out,'/outs/fragments.tsv.gz')
fragments      <- CreateFragmentObject(path = fragments.path,
                                       cells = metadata$barcode,
                                       verbose = TRUE,validate.fragments = TRUE)


############################ Create peak counts matrix
cat("*** Creating bins/peaks matrices \n")

peaks_file = paste0('results/',args$sample,'/macs/broad/',args$sample,'_peaks.broadPeak')
if (!file.exists(peaks_file)) {stop(paste0("Peaks file does not exist:: ", peaks_file))}
peaks <- rtracklayer::import(peaks_file)


counts.matrix.peaks <- FeatureMatrix(fragments = fragments,
                                     features = peaks,
                                     cells = metadata$barcode)


counts.matrix.bins <- GenomeBinMatrix(fragments = fragments,
                                      genome = setNames(getChromInfoFromUCSC(genome_version)[,2],getChromInfoFromUCSC(genome_version)[,1]),
                                      binsize = args$window,
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
# seurat_object            <- RenameCells(object = seurat_object,add.cell.id = args$sample)

new.metadata <- unlist(config$samples[[args$sample]])
for (x in seq(new.metadata)) {
  seurat_object <- AddMetaData(object = seurat_object,metadata = new.metadata[x],col.name = names(new.metadata)[x])
}

########## Create Gene activity matrix
genebodyandpromoter.coords.flat <- load_ensembl_annot(genome_version)

gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat ,
                                 cells = colnames(seurat_object))

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)
rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]

gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]


seurat_object[['GA']] <- CreateAssayObject(counts = gene.matrix)

######### Save the object
cat("*** Save the object \n")
dir.create(args$out_prefix,recursive = TRUE)
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

########################## Try Clustering using Seurat
cat("*** Clustering and dimensionality reduction \n")

DefaultAssay(seurat_object) <- paste0('bins_',args$window)
# DefaultAssay(seurat_object) <- "peaks"

# The dimensionality reduction might fail, especially for low complexity datasets, so put the code into try to still save the objects
try({
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object,min.cutoff = 'q0')
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
  #algorithm = "leiden",
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

