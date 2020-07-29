# Preprocessing of the H3K4me3 dataset


print(.libPaths())

if (!'Signac' %in% rownames(installed.packages())){
  cat("Installing missing package Signac \n")
  install.packages('Signac',repos='http://cran.us.r-project.org')
}

if (!'harmony' %in% rownames(installed.packages())){
  cat("Installing missing package harmony \n")
  devtools::install_github("immunogenomics/harmony")
}



cat("*** Loading libraries")
suppressMessages({
library(argparse);
library(Seurat);
library(Signac);
library(EnsDb.Mmusculus.v79);
library(viridis);
library(ggplot2);
library(dplyr);
library(BSgenome.Mmusculus.UCSC.mm10)
# library(BSgenome.Mmusculus.UCSC.mm10);
  }
)

set.seed(1234)

########### Arguments parser

parser <- ArgumentParser()

parser$add_argument("-s", "--sample", type="character", default='foo', 
                    help="sample name [as in config file key]")

parser$add_argument("-c", "--config", type="character", default='config/config.yaml', 
                    help="maximum number of reads in cell")

parser$add_argument("-o", "--out_prefix", type="character", default="100000_UMI", 
                    help="folder for the output in clustering_snakemake folder")

parser$add_argument("-w", "--window", type="integer", default=10000, 
                    help="width of the genome window (if -t bins)")

args <- parser$parse_args()

# args <- list()
# args$sample     <- "H3K27me3_N1"
# args$config     <- "scCut-Tag_2020/config/config.yaml"
# args$out_prefix <- "results/H3K27me3_N1/cell_picking/5000/"
# args$window     <- 5000


print(args)

config <- yaml::read_yaml(args$config)


cutoff_reads_min            = config$samples[[args$sample]]$clustering_params$min_reads_log10
cutoff_reads_max            = config$samples[[args$sample]]$clustering_params$max_reads_log10
cutoff_peak_percentage_low  = config$samples[[args$sample]]$clustering_params$min_peaks_ratio
cutoff_peak_percentage_high = config$samples[[args$sample]]$clustering_params$max_peaks_ratio

ndim = 30
window  = args$window

fragments <- paste0(config$samples[[args$sample]]$cellranger_out,'/outs/fragments.tsv.gz')
assay = "peaksMB"


if (!file.exists(fragments)) {stop(paste0("Fragments file does not exist: ",fragments))}

#### Read annotation

cat("*** Loading mm10 annotation \n")
ensdb = EnsDb.Mmusculus.v79

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



########## Filter the barcodes
cat("*** Reading barcode statistics files \n")

all_barcodes_file  <- paste0('results/',args$sample,'/barcode_statistics/all_barcodes.txt')
peak_barcodes_file <- paste0('results/',args$sample,'/barcode_statistics/peaks_barcodes_broad.txt')
metadata_file      <- paste0(config$samples[[args$sample]]$cellranger_out,'/outs/singlecell.csv')

metadata = read.csv(metadata_file, header = 1)
metadata = metadata[2:nrow(metadata),]
metadata$logUMI = log10(metadata$passed_filters + 1)
metadata$promoter_ratio = (metadata$promoter_region_fragments+1) / (metadata$passed_filters + 1)
metadata$peak_region_ratio = (metadata$peak_region_fragments+1) / (metadata$passed_filters + 1)

all_barcodes <- read.table(file=all_barcodes_file)
peak_barcodes <- read.table(file=peak_barcodes_file)
bcd <- merge(all_barcodes,peak_barcodes,by="V2")  
colnames(bcd) <- c("barcode","all_unique_MB","peak_MB")
bcd$peak_ratio_MB <- bcd$peak_MB/bcd$all_unique_MB
bcd$sample <- args$sample

# Merge 10x metadata with barcode statistics
metadata <- merge(metadata,bcd,by='barcode')


metadata$is__cell_barcode <- as.factor(metadata$is__cell_barcode)

################ MB filtering
cat("*** Filtering cells \n")
metadata[,"passedMB"] <- FALSE
metadata[metadata$all_unique_MB > 10^cutoff_reads_min &
         metadata$all_unique_MB < 10^cutoff_reads_max &
         metadata$peak_ratio_MB > cutoff_peak_percentage_low &
         metadata$peak_ratio_MB < cutoff_peak_percentage_high,"passedMB"] <- TRUE


################ Cell picking scatterplot nreads ~ percent in peaks

p1 <- ggplot(data = metadata,aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=is__cell_barcode)) + 
  geom_point(size=0.1) + 
  scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(as.numeric(as.character(metadata$is__cell_barcode)))),"FALSE")) +
  theme(legend.position="bottom",text=element_text(size=26)) 

p2 <- ggplot(data = metadata,aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=is__cell_barcode)) + 
  geom_point(size=0.1) +
  scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(as.numeric(as.character(metadata$is__cell_barcode)))),NA)) + 
  theme(legend.position="bottom",text=element_text(size=26))

ggsave(plot = p1+p2,
       filename=paste0(args$out_prefix,'cells_10x.png'),width = 20,height = 10,units = 'in')


p1 <- ggplot(data = metadata) + 
      geom_point(aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=passedMB),size=0.1) + 
      geom_hline(yintercept = c(cutoff_peak_percentage_high,cutoff_peak_percentage_low)) + 
      geom_vline(xintercept = c(cutoff_reads_min,cutoff_reads_max)) + 
      scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(metadata$passedMB)),"FALSE")) +
      theme(legend.position="bottom",text=element_text(size=26))
  
p2 <- ggplot(data = metadata) + 
      geom_point(aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=passedMB),size=0.1) +
      scale_color_manual(values=c("black","gold")) + 
      theme(legend.position="bottom",text=element_text(size=26))

ggsave(plot = p1+p2,
       filename=paste0(args$out_prefix,'cells_picked.png'),width = 20,height = 10,units = 'in')


################# Export bw selected / unselected
cat("*** Reading fragments file \n")

fragments_gr      <- rtracklayer::import(fragments,format = "bed")


cat("*** Exorting merged bw files \n")
barcode_pass   <- metadata$barcode[metadata$passedMB]
barcode_nopass <- metadata$barcode[!metadata$passedMB]
  
fragments.pass   <- fragments_gr[fragments_gr$name %in% barcode_pass]
fragments.nopass <- fragments_gr[fragments_gr$name %in% barcode_nopass]
  
cat("*** Calculating coverage \n")
coverage.pass   <- GenomicRanges::coverage(fragments.pass)
coverage.nopass <- GenomicRanges::coverage(fragments.nopass)
  
cat("*** Normalizing \n")
coverage.pass <- coverage.pass/length(fragments.pass)
coverage.nopass <- coverage.nopass/length(fragments.nopass)
  
cat("*** Exporting \n")
rtracklayer::export(object=coverage.pass,  con = paste0(args$out_prefix,'cells_picked.bw'))
rtracklayer::export(object=coverage.nopass,con = paste0(args$out_prefix,'cells_not_picked.bw'))

############################ Filter the dataset

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

