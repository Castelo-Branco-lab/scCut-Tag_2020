#### Read annotation

load_ensembl_annot <- function() {
  library(EnsDb.Mmusculus.v79)
  library(ensembldb)
  
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
  return(genebodyandpromoter.coords.flat)
}


get_path_to_script <- function(){
  args <- commandArgs(trailingOnly=FALSE)
  x    <- grep(pattern="--file=",x=args)
  path <- args[x]
  path <- gsub("--file=","",path)
  return(path)
}