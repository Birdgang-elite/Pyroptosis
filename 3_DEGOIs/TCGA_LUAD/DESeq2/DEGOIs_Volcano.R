rm(list=ls())
setwd("C:\\Users\\xiaog\\Desktop\\Pyroptosis\\3_DEGOIs\\TCGA_LUAD\\DESeq2")

DEGsup <- read.table("DEGsup_DESeq2.txt",sep="\t",header=T,check.names=F) 
DEGsdown <- read.table("DEGsdown_DESeq2.txt",sep="\t",header=T,check.names=F) 
GOIs <- read.table("GOIs.txt",sep="\t",header=T,check.names=F) 
colnames(GOIs) <- "gene"

DEGOIsup <- merge(DEGsup,GOIs, by="gene")
write.table(DEGOIsup, file="TCGA_DESeq2_DEGOIsup.txt",sep="\t",quote=F,row.names=F)
DEGOIsdown <- merge(DEGsdown,GOIs, by="gene")
write.table(DEGOIsdown, file="TCGA_DESeq2_DEGOIsdown.txt",sep="\t",quote=F,row.names=F)



#Volcano
foldChange =1
padj = 0.05

DEGs <- read.table("DEGs_DESeq2.txt",sep="\t",header=T,check.names=F)
DEGOIs <- merge(DEGs,GOIs, by="gene")

#Volcano plot
pdf(file="UCSC_vol.pdf")
xMax=max(-log10(DEGOIs$padj))+1
yMax=5
plot(-log10(DEGOIs$padj), DEGOIs$log2FoldChange, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=DEGOIs[DEGOIs$padj<padj & DEGOIs$log2FoldChange>foldChange,]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=DEGOIs[DEGOIs$padj<padj & DEGOIs$log2FoldChange<(-foldChange),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

