library(Seurat)
library(SeuratDisk)

#args <- "/data/proj/GCB_MB/CT/git_test/results/H3K4me3/clustering/01.clustering.Rds"

args <- commandArgs(trailingOnly=TRUE)
prefix   <- dirname(args[1])
antibody <- unlist(strsplit(prefix,split = "/"))[8]


# Load Seurat object
object <- readRDS(args[1])

# Create output h5 filename
out_file <- paste(prefix,paste0(antibody,"_matrix_h5"),sep="/")

# Export all assays
SaveH5Seurat(object=object,filename=out_file)



