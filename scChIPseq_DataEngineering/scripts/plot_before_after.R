library(ggplot2)
library(dplyr)

odir = as.character(commandArgs(TRUE)[1])
name = as.character(commandArgs(TRUE)[2])

before=read.table(sep = "\t", paste0(odir,"/",name,"_flagged.bc_file"),header =  T)

after=read.table(sep = "\t", paste0(odir,"/",name,"_flagged_rmDup.bc_file"),header =  T)

reads_before_after = dplyr::right_join(before,after,by="value")

colnames(reads_before_after) = c("Barcode","ReadsBefore_rmDup","ReadsAfter_rmDup")
reads_before_after$ReadsBefore_rmDup = as.numeric(reads_before_after$ReadsBefore_rmDup)
reads_before_after$ReadsAfter_rmDup = as.numeric(reads_before_after$ReadsAfter_rmDup)

lm = lm(reads_before_after$ReadsAfter_rmDup ~ reads_before_after$ReadsBefore_rmDup)

pdf(paste0(as.character(odir),'/before_vs_after_rmDup.pdf'))
corr = round(cor(reads_before_after$ReadsBefore_rmDup, reads_before_after$ReadsAfter_rmDup, method="pearson"), 4)

p <- ggplot(reads_before_after, aes(x=ReadsAfter_rmDup, y=ReadsBefore_rmDup)) +
  geom_point(alpha=0.5, colour="steelblue") +  scale_x_log10() + scale_y_log10() +
  labs(title=paste0("Total counts /cell before vs after ( #readsAfter =  ", round(lm$coefficients[2],3), " x #readsBefore + ",round(lm$coefficients[1],1), " - R² = ",round(summary(lm)$r.squared,3),"- corr = ",corr,")"), x="#reads after", y="#reads before") +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed")  + theme(plot.title = element_text(size=10))
plot(p)
dev.off()

