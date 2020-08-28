#!/usr/bin/env R

# Usage 
# File1 column File2 column File3 column ... FileX column out_file
# Column is 1-indexed

args     <- commandArgs(trailingOnly=TRUE)

out_file <- args[length(args)]
args     <- args[-length(args)]

files  <- args[seq(1,length(args),2)]
column <- args[seq(1,length(args),2)+1]
column <- as.numeric(as.character(column))

dfs <- lapply(seq(files),function(x){
  df <- read.table(file = files[x])
  colnames(df)[column[x]] <- 'merge_by'
  df
})

df_final <- Reduce(function(df1, df2) merge(df1, df2, by = "merge_by"),dfs)
df_final <- df_final[order(as.numeric(gsub('loop_','',df_final$merge_by))),]
df_final <- df_final[,-1]

write.table(x = df_final,file = out_file,quote = FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
