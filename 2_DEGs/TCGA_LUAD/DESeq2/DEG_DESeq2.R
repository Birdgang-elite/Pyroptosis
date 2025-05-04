rm(list = ls())
setwd("C:\\Users\\xiaog\\Desktop\\LUAD\\2_DEGs_of_each_cohort\\TCGA_LUAD\\DESeq2") 
library("limma")
library("edgeR")
library("DESeq2")

##Set abs(logFC) and Padj
log2FC=1
Padj=0.05

##Build the expression matrix and grouplist.
exprs <- read.table("mRNA.symbol.txt",sep="\t",header=T,check.names=F) #59 Normal in front of 535 LUAD
exprs1 <- avereps(exprs[,-1],ID = exprs$id) #对重复基因名取平均表达量，然后将基因名作为行名
exprs2 <- exprs1[rowMeans(exprs1)>1,] #去除低表达的基因
LUAD <- read.table("LUAD.txt",sep="\t",header=T,check.names=F)
Normal <- read.table("Normal.txt",sep="\t",header=T,check.names=F)
grouplist <- factor(c(rep("Normal",nrow(Normal)),rep("LUAD",nrow(LUAD))),levels = c('Normal','LUAD'))

##DESeq2
dat <- exprs2
dat <- apply(dat,2,as.integer)#整数转换
rownames(dat) <- rownames(exprs2)
colData <- data.frame(row.names =colnames(dat), 
                      condition=grouplist)
dds <- DESeqDataSetFromMatrix(
  countData = dat,
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition",rev(levels(grouplist))))
resOrdered <- res[order(res$pvalue),] # 按照P值排序
DEGs <- as.data.frame(resOrdered)
DEGs = na.omit(DEGs)
write.table(DEGs, file="DEGs_DESeq2.txt",sep="\t",quote=F,row.names=T)
DEGsup=DEGs[DEGs$padj<Padj & DEGs$log2FoldChange>log2FC,]
write.table(DEGsup, file="DEGsup_DESeq2.txt",sep="\t",quote=F,row.names=T)
DEGsdown=DEGs[DEGs$padj<Padj & DEGs$log2FoldChange<(-log2FC),] 
write.table(DEGsdown, file="DEGsdown_DESeq2.txt",sep="\t",quote=F,row.names=T)

