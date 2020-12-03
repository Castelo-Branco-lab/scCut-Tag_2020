#!/usr/bin/env R

library(reshape2)
library(ggplot2)
library(ggthemes)

#args      <- commandArgs(trailingOnly=TRUE)
args = c('results/other_datasets/fingerprint_analysis/scCT_fingerprint.txt',
         'results/other_datasets/fingerprint_analysis/Kaya_okur_fingerprint.txt',
         'results/other_datasets/fingerprint_analysis/Grosselin_fingerprint.txt')


exp.names <- unlist(lapply(strsplit(basename(args),'_'),'[',1))


to.plot <- c(c("K562_H3K27me3_iCell8.bam","K562_H3K4me2_iCell8.bam"),
             c("Astrocytes","Microglia","mOL","Neurons","Neurons","OEC","OPC","VLMC","/outs/possorted_bam.bam"),
             c("possorted_SRR7536860.bam"))

scCT_bulk <- "scCT_all_clusters.bam"

# Read files
df.ls <-lapply(args,function(x){
  NAMES       <- read.table(file=x,nrow=1,stringsAsFactors=F)
  D           <- as.data.frame(read.table(file=x,skip=2,stringsAsFactors=F))
  colnames(D) <- NAMES
  D[,grep(paste(to.plot,collapse="|"),colnames(D)),drop=FALSE]
})


# Make cumsum
df.ls.cumsum <- lapply(df.ls,function(x){
  x <- apply(x,2,sort)
  x <- apply(x,2,function(y){
    y <- cumsum(y)/sum(y)
    return(y)
  })
  x
})

colnames(df.ls.cumsum[[1]]) <- c('scCT_bulk_N1','scCT_bulk_N2','scCT_bulk_N3','scCT_bulk_N4',colnames(df.ls.cumsum[[1]])[5:length(colnames(df.ls.cumsum[[1]]))])


# Melt
df.ls.cumsum.melt <- lapply(df.ls.cumsum,function(x){
  x <- melt(x)
  x
})

# add names
df.ls.cumsum.melt <- lapply(1:length(df.ls.cumsum),function(x){
  df.ls.cumsum.melt[[x]]$group <- exp.names[x]
  df.ls.cumsum.melt[[x]]
})

df      <- do.call('rbind',df.ls.cumsum.melt)

# Fix names for plotting
df[grep('scCT_bulk_N',df$Var2),'group']   <- 'scCT_merged_replicate'
df[df$group=='scCT','group']               <- 'scCT_single_cluster'

#antibody column
df$antibody <- "H3K27me3"
df[grep("H3K4me2",df$Var2),'antibody'] <- "H3K4me2"


df$Var2 <- basename(as.character(df$Var2))
df$Var2 <- gsub(".bam","",df$Var2)


ggplot(data=df) + geom_line(aes(x=Var1,y=value,col=group,group=Var2,lty=antibody)) + 
  coord_cartesian(xlim=c(200000,max(df$Var1))) + 
  theme_few()
