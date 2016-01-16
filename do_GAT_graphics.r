library(ggplot2)
myData<-read.table("GAT_GM_consensus.tsv", header=T)
ggplot(myData, aes(y=l2fold, x=reorder(annotation,l2fold)), fill=annotation) + geom_bar(stat='identity',colour="black",fill="lightgreen") + coord_flip() + ggtitle("GAT (consensus peakset)") + ylab("log2(fold change)") + xlab("Annotation (Ensembl 67)") + theme(plot.title = element_text(face = "bold", size=15)) + theme(axis.text.y = element_text(family = "sans", face = "bold", size = 12)) + theme(axis.text.x = element_text(family = "sans", face = "bold", size = 12)) -> p
pdf('test.pdf')
print(p)
dev.off()



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
