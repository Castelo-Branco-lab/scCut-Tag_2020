
clusterSeurat <- function(brain,ndim) {
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
  
  p1 <- DimPlot(brain,group.by = 'ident',pt.size = 0.2) + NoLegend() + ggtitle(paste0("ndim= ",ndim))
  return(p1)
}

clusterMonocle <- function(brain,assay,ndim){
  cds <- new_cell_data_set(expression_data = brain[[assay]]@counts, cell_metadata = brain@meta.data)
  cds <- preprocess_cds(cds, num_dim = ndim, method = "LSI")
  cds <- reduce_dimension(cds,preprocess_method = "LSI", reduction_method="UMAP")
  cds <- cluster_cells(cds)
  cds <- invisible(learn_graph(cds))
  
  p1 <- plot_cells(cds)
  p1  
}



#plots <- lapply(c(3:15,20,25,30),function(x){
#  clusterSeurat(brain,x)
#  })

#g <- do.call("grid.arrange",plots)
#ggsave(filename = "01.Seurat_peaks_ndim_matrix.pdf",plot = g,width=12,heigh=12)

  

