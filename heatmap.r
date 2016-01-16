library(RColorBrewer)
library(gplots)
x=read.table("matrix.dat", header=TRUE)
mat=data.matrix(x)

pdf("heatmap.pdf", height=10, width=10)
heatmap.2(mat,
Rowv=TRUE,
Colv=TRUE,
#    dendrogram= c("none"),
distfun = dist,
hclustfun = hclust,
xlab = "X data", ylab = "Y data",
key=TRUE,
keysize=1,
trace="none",
density.info=c("none"),
margins=c(10, 8),
col=brewer.pal(10,"PiYG")
#    col=redgreen(75),
)
dev.off()
