#!/usr/bin/env Rscript 

args <- commandArgs(TRUE);
# args <- c("/groups/berger/user/pierre.bourguet/genomics/Araport11/tmp2")
# args <- c("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3_polyA/DEGs_cmt23_vs_WT/")

x <- read.delim(args[1], header=F, sep='\t', quote="", dec=".", comment.char="")
names(x) <- c("chr", "start", "end", "nb_DEG")
x$chr <- substr(x$chr,4,4)
# plot log2FC mean for 100kb windows along chromosomes
pdf(paste0(dirname(args[1]), "/100kb_plot_", substr(basename(args[1]), 1, nchar(basename(args[1]))-21) ,".pdf"), width = 12)
par(mfrow=c(2,3))
for (k in 1:5) {
  par(pin = c(max(x$end[x$chr == k], na.rm = T) / 10000000, 1.899979))
  plot(x$end[x$chr == k], x$nb_DEG[x$chr == k]
       , ylim = c(floor(min(x$nb_DEG, na.rm = T)),ceiling(max(x$nb_DEG, na.rm = T))) #the max value to be plotted is extracted, ylim is defined as the next integer to this max
       , type = "h"
       , ylab = "nb of DEGs per 100kb window"
       , xlab = NA
       , xaxt = "n"
       , mgp = c(2.5,1,0) #margin line for axis title, axis labels and axis line
       , mar = c(1,0,2,0)
       , cex.axis = 1
       , cex.lab = 1.2
  )
  mtext(paste("chromosome", k), side = 1 #add "chromosome X" under x axis
        , cex = 0.8, line = 1 #magnification & position of the text
  )
  axis(side = 1, at = c(x$start[x$chr == k][1], x$end[x$chr == k][nrow(x[x$chr == k,])]  ), labels = c(x$start[x$chr == k][1], x$end[x$chr == k][nrow(x[x$chr == k,])]  ), cex.axis = 0.9)
}
title(main = substr(basename(args[1]), 1, nchar(basename(args[1]))-21), outer = TRUE, line = -3, cex.main = 2)

graphics.off()
