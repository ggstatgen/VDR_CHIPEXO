library(gplots)

#test
vdr <- read.table("/net/isi-scratch/giuseppe/VDR/RAW/34_SVR11_hg19_LOW_QUALITY/02_BAM_BOWTIE/d_correlation_all_q0.01/d_loglr/wc_R.dat",header=TRUE, row.names=1)
pdf("/net/isi-scratch/giuseppe/VDR/RAW/34_SVR11_hg19_LOW_QUALITY/02_BAM_BOWTIE/d_correlation_all_q0.01/d_loglr/heatmap.pdf", height=10, width=10)
heatmap.2(as.matrix(vdr),trace="none")
dev.off()
