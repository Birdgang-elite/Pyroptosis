rm(list = ls())
path <- getwd()
setwd(path)

library(ggpubr)
library(tidyverse)
library(reshape2)

senstivity <- read.csv('DrugPredictions.csv',
                       header = T, sep = ',', check.names = F,row.names = 1)
colnames(senstivity) <- gsub("(.*)\\_(\\d+)","\\1",colnames(senstivity))

senstivity$V199 <- substr(row.names(senstivity),1,15)
senstivity1 <- senstivity[!duplicated(senstivity$V199),]
rownames(senstivity1) <- senstivity1[,dim(senstivity)[2]]
senstivity <- as.data.frame(senstivity1[,-199])
rownames(senstivity) <- rownames(senstivity1)

#读入风险文件
risk <- read.table('UCSC_LUAD_exprs_RS_all.txt', header = T, sep = "\t", check.names = F, row.names = 1)
#合并
sameSample <- intersect(row.names(risk),row.names(senstivity))
risk <- risk[sameSample,"Risk_score_LH",drop=F]
senstivity <- senstivity[sameSample,,drop=F]
senstivity[is.na(senstivity)] <- 0
rt <- cbind(risk,senstivity)

#设置比较组
rt$Risk_score_LH <- factor(rt$Risk_score_LH, levels=c(0,1),labels=c('Risk-low','Risk-high'))
type <- levels(factor(rt[,"Risk_score_LH"]))
comp <- combn(type,2)
my_comparisions <- list()
for(i in 1:ncol(comp)){my_comparisions[[i]] <- comp[,i]}


#提取显著差异的药物
sigGene <- c()
for(i in colnames(rt)[2:(ncol(rt))]){
  if(sd(rt[,i])<0.05){next}
  wilcoxTest = wilcox.test(rt[,i]  ~ rt[,"Risk_score_LH"])
  pvalue = wilcoxTest$p.value
  if(wilcoxTest$p.value < 0.001){                 ######根据需要设定p.value
    sigGene = c(sigGene,i)
  }
}
sigGene <- c(sigGene,"Risk_score_LH")
rt <- rt[,sigGene]

#把数据转换成ggplot2输入文件
rt <- melt(rt,id.vars <- c("Risk_score_LH"))
colnames(rt) <- c("Risk_score_LH", "Gene", "Expression")

#设置比较组
group <- levels(factor(rt$Risk_score_LH))

comp <- combn(group,2)
my_comparisions <- list()
for(j in 1:ncol(comp)){my_comparisions[[j]] <- comp[,j]}

#绘制箱线图
boxplot <- ggboxplot(rt, x="Gene", y="Expression", fill = "Risk_score_LH",
                     xlab = "",
                     ylab = "Imputed senstivity score",
                     legend.title = "Risk_score_LH",
                     width = 0.8,
                     palette = c('DodgerBlue1','Firebrick2'))+
  rotate_x_text(50)+
  stat_compare_means(aes(group = Risk_score_LH),
                     method = 'wilcox.test',
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols = c("***","**","*","ns")),label = "p.signif")+
  theme(axis.text=element_text(face = "bold.italic",colour="#441718",size=16),
        axis.title =element_text(face ="bold.italic",colour='#441718',size=16),
        axis.line =element_blank(),
        plot.title=element_text(face ="bold.italic",colour ='#441718',size=16),
        legend.text=element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A790",size=1.5, linetype ="solid"),
        panel.background=element_rect(fill ="#F1F6FC"),
        panel.grid.major=element_line(color = "#CFD3D6",size =.5,linetype ="dotdash"),
        legend,title=element_text(face ="bold.italic",size =13)
  )

#输出图片
pdf(file = "drugSenstivity.pdf",width = 20, height = 8)
print(boxplot)
dev.off()


#分别绘制sigGene中药物敏感性的柱状图
for(i in 1:length(sigGene)-1){
  p <- ggplot(subset(rt,Gene == sigGene[i]),aes(Gene,Expression, fill = Risk_score_LH)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    #scale_fill_manual(values = palette1[c(2,4)])+ 
    theme_bw() + 
    labs(x = NULL, y = "Imputed senstivity score") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=0,hjust = 1),
          axis.text = element_text(color = "black",size = 12))+
    stat_compare_means()
  ggsave(paste0(sigGene[i],'-Bar.pdf'),plot = p,width = 7,height = 7) 
}



####DEGs
rm(list = ls())
path <- getwd()
setwd(path)

library(ggpubr)
library(tidyverse)
library(reshape2)
library(limma)

senstivity <- read.csv('DrugPredictions.csv',
                       header = T, sep = ',', check.names = F,row.names = 1)
colnames(senstivity) <- gsub("(.*)\\_(\\d+)","\\1",colnames(senstivity))

senstivity$V199 <- substr(row.names(senstivity),1,15)
senstivity1 <- senstivity[!duplicated(senstivity$V199),]
rownames(senstivity1) <- senstivity1[,dim(senstivity)[2]]
senstivity <- as.data.frame(senstivity1[,-199])
rownames(senstivity) <- rownames(senstivity1)
GDSC2 <- rownames_to_column(senstivity,'Sample')

UCSC_LUAD_exprs_RS_all <- "UCSC_LUAD_exprs_RS_all.txt" %>%
  read.table(sep = "\t",header=T,check.names=F)

UCSC_LUAD_GDSC2 <- UCSC_LUAD_exprs_RS_all %>%
  subset(select = c(Sample,Risk_score,Risk_score_LH)) %>%
  merge(GDSC2,by='Sample')


RS_High_GDSC2 <- subset(UCSC_LUAD_GDSC2,Risk_score_LH == 1) %>%
  t() %>% as.data.frame()
colnames(RS_High_GDSC2) <- RS_High_GDSC2[1,]
RS_High_GDSC2 <- RS_High_GDSC2[-1,]
RS_High_GDSC2 <- rownames_to_column(RS_High_GDSC2,'Drug')

RS_Low_GDSC2 <- subset(UCSC_LUAD_GDSC2,Risk_score_LH == 0) %>%
  t() %>% as.data.frame()
colnames(RS_Low_GDSC2) <- RS_Low_GDSC2[1,]
RS_Low_GDSC2 <- RS_Low_GDSC2[-1,]
RS_Low_GDSC2 <- rownames_to_column(RS_Low_GDSC2,'Drug')

All_GDSC2 <- merge(RS_High_GDSC2,RS_Low_GDSC2,by = 'Drug')  
rownames(All_GDSC2) <- All_GDSC2[,1]
All_GDSC2 <- All_GDSC2[,-1]
All_GDSC2_2 <- as.data.frame(lapply(All_GDSC2,as.numeric))
rownames(All_GDSC2_2) <- rownames(All_GDSC2)
colnames(All_GDSC2_2) <- colnames(All_GDSC2)

grouplist <- factor(c(rep("RS_H",ncol(RS_High_GDSC2)-1),rep("RS_L",ncol(RS_Low_GDSC2)-1)),levels = c('RS_H','RS_L'))
design <- model.matrix(~0+grouplist)
colnames(design) <- c("RS_H","RS_L")
fit <- lmFit(All_GDSC2_2, design)
cont.matrix <- makeContrasts(RS_HvsRS_L=RS_H-RS_L, levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
DEGs <- topTable(fit2, n = Inf) %>%
  rownames_to_column('Factor')
write.table(DEGs, file="DEGs_GDSC2.txt",sep="\t",quote=F,row.names=F)

