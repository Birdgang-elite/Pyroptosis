rm(list=ls())
path <- getwd()
setwd(path) 
library(limma)
library(ggpubr)
library(tidyverse)


#读取TIDE数据
tide <- read.csv("TCGA_LUAD_TIDE.csv",header = T, sep = ',', check.names = F)
tide <- tide %>%
  mutate(Sample = substr(tide$Patient,1,15))
tide <- subset(tide,select = c('Sample','TIDE'))

UCSC_LUAD_exprs_RS_all <- read.table("UCSC_LUAD_exprs_RS_all.txt",header = T, sep = "\t", check.names = F)

UCSC_LUAD_RS_TIDE <- merge(UCSC_LUAD_exprs_RS_all,tide,by='Sample') 
UCSC_LUAD_RS_TIDE$Risk_score_LH <- factor(UCSC_LUAD_RS_TIDE$Risk_score_LH,levels=c(0,1),labels=c('Risk-low','Risk-high'))


#柱状图
p <- ggplot(UCSC_LUAD_RS_TIDE,aes(Risk_score_LH,TIDE, fill = Risk_score_LH)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  #scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "TIDE") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=0,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means()
ggsave('TIDE-Bar.pdf',plot = p,width = 7,height = 7)  


#小提琴图
p1 <- ggviolin(UCSC_LUAD_RS_TIDE,x='Risk_score_LH',y='TIDE',fill = 'Risk_score_LH',
               xlab = '',ylab = 'TIDE',
               palette = c('#0066FF','#FF0000'),
               legend.title = 'Risk_score_LH',
               add = 'boxplot',add.params = list(fill='white'))+
  stat_compare_means()
ggsave('TIDE-violin.pdf',plot = p1,width = 7,height = 7)  


#Pearson相关性分析的散点图
p2 <- ggplot(UCSC_LUAD_RS_TIDE,aes(x = Risk_score, y = TIDE))+
  geom_point()+
  geom_smooth(method = "lm",color = "black", fill = "lightgray")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x = "Risk score", y = "TIDE") +
  stat_cor(method = "pearson",
           label.x = 1.5,label.y = 2.5)
ggsave('TIDE_correlation.pdf',plot = p2,width = 7,height = 7)  