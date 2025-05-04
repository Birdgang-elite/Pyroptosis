################################################
##请选取所需的数据集，无需调整本代码。##########
################################################
rm(list=ls())
path <- getwd()
setwd(path)
library(tidyverse)


#GSE3141
GSE3141_LUAD_exprs <- paste0(path,'\\GSE3141\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
GSE3141_LUAD_exprs=data.frame(t(GSE3141_LUAD_exprs))
colnames(GSE3141_LUAD_exprs) <- GSE3141_LUAD_exprs[1,]
GSE3141_LUAD_exprs <- GSE3141_LUAD_exprs[-1,]
GSE3141_LUAD_exprs$Sample <- rownames(GSE3141_LUAD_exprs)

GSE3141_LUAD <- paste0(path,'\\GSE3141\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(GSE3141_LUAD) <- GSE3141_LUAD[,1]
GSE3141_LUAD_os <- GSE3141_LUAD[,c('Sample','os_time','os_status')]

GSE3141_LUAD_exprs_os <- merge(GSE3141_LUAD_os,GSE3141_LUAD_exprs,by = 'Sample')
GSE3141_LUAD_exprs_os1 <- GSE3141_LUAD_exprs_os[GSE3141_LUAD_exprs_os$os_time>30,]
GSE3141_LUAD_exprs_os2 <- na.omit(GSE3141_LUAD_exprs_os1)
colnames(GSE3141_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(GSE3141_LUAD_exprs_os2, file=paste0(path,'\\GSE3141\\GSE3141_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)

#GSE30219
GSE30219_LUAD_exprs <- paste0(path,'\\GSE30219\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
GSE30219_LUAD_exprs=data.frame(t(GSE30219_LUAD_exprs))
colnames(GSE30219_LUAD_exprs) <- GSE30219_LUAD_exprs[1,]
GSE30219_LUAD_exprs <- GSE30219_LUAD_exprs[-1,]
GSE30219_LUAD_exprs$Sample <- rownames(GSE30219_LUAD_exprs)

GSE30219_LUAD <- paste0(path,'\\GSE30219\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(GSE30219_LUAD) <- GSE30219_LUAD[,1]
GSE30219_LUAD_os <- GSE30219_LUAD[,c('Sample','os_time','os_status')]

GSE30219_LUAD_exprs_os <- merge(GSE30219_LUAD_os,GSE30219_LUAD_exprs,by = 'Sample')
GSE30219_LUAD_exprs_os1 <- GSE30219_LUAD_exprs_os[GSE30219_LUAD_exprs_os$os_time>30,]
GSE30219_LUAD_exprs_os2 <- na.omit(GSE30219_LUAD_exprs_os1)
colnames(GSE30219_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(GSE30219_LUAD_exprs_os2, file=paste0(path,'\\GSE30219\\GSE30219_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)

#GSE31210
GSE31210_LUAD_exprs <- paste0(path,'\\GSE31210\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
GSE31210_LUAD_exprs=data.frame(t(GSE31210_LUAD_exprs))
colnames(GSE31210_LUAD_exprs) <- GSE31210_LUAD_exprs[1,]
GSE31210_LUAD_exprs <- GSE31210_LUAD_exprs[-1,]
GSE31210_LUAD_exprs$Sample <- rownames(GSE31210_LUAD_exprs)

GSE31210_LUAD <- paste0(path,'\\GSE31210\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(GSE31210_LUAD) <- GSE31210_LUAD[,1]
GSE31210_LUAD_os <- GSE31210_LUAD[,c('Sample','os_time','os_status')]

GSE31210_LUAD_exprs_os <- merge(GSE31210_LUAD_os,GSE31210_LUAD_exprs,by = 'Sample')
GSE31210_LUAD_exprs_os1 <- GSE31210_LUAD_exprs_os[GSE31210_LUAD_exprs_os$os_time>30,]
GSE31210_LUAD_exprs_os2 <- na.omit(GSE31210_LUAD_exprs_os1)
colnames(GSE31210_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(GSE31210_LUAD_exprs_os2, file=paste0(path,'\\GSE31210\\GSE31210_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)

#GSE41271
GSE41271_LUAD_exprs <- paste0(path,'\\GSE41271\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
GSE41271_LUAD_exprs=data.frame(t(GSE41271_LUAD_exprs))
colnames(GSE41271_LUAD_exprs) <- GSE41271_LUAD_exprs[1,]
GSE41271_LUAD_exprs <- GSE41271_LUAD_exprs[-1,]
GSE41271_LUAD_exprs$Sample <- rownames(GSE41271_LUAD_exprs)

GSE41271_LUAD <- paste0(path,'\\GSE41271\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(GSE41271_LUAD) <- GSE41271_LUAD[,1]
GSE41271_LUAD_os <- GSE41271_LUAD[,c('Sample','os_time','os_status')]

GSE41271_LUAD_exprs_os <- merge(GSE41271_LUAD_os,GSE41271_LUAD_exprs,by = 'Sample')
GSE41271_LUAD_exprs_os1 <- GSE41271_LUAD_exprs_os[GSE41271_LUAD_exprs_os$os_time>30,]
GSE41271_LUAD_exprs_os2 <- na.omit(GSE41271_LUAD_exprs_os1)
colnames(GSE41271_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(GSE41271_LUAD_exprs_os2, file=paste0(path,'\\GSE41271\\GSE41271_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)

#GSE50081
GSE50081_LUAD_exprs <- paste0(path,'\\GSE50081\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
GSE50081_LUAD_exprs=data.frame(t(GSE50081_LUAD_exprs))
colnames(GSE50081_LUAD_exprs) <- GSE50081_LUAD_exprs[1,]
GSE50081_LUAD_exprs <- GSE50081_LUAD_exprs[-1,]
GSE50081_LUAD_exprs$Sample <- rownames(GSE50081_LUAD_exprs)

GSE50081_LUAD <- paste0(path,'\\GSE50081\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(GSE50081_LUAD) <- GSE50081_LUAD[,1]
GSE50081_LUAD_os <- GSE50081_LUAD[,c('Sample','os_time','os_status')]

GSE50081_LUAD_exprs_os <- merge(GSE50081_LUAD_os,GSE50081_LUAD_exprs,by = 'Sample')
GSE50081_LUAD_exprs_os1 <- GSE50081_LUAD_exprs_os[GSE50081_LUAD_exprs_os$os_time>30,]
GSE50081_LUAD_exprs_os2 <- na.omit(GSE50081_LUAD_exprs_os1)
colnames(GSE50081_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(GSE50081_LUAD_exprs_os2, file=paste0(path,'\\GSE50081\\GSE50081_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)

#UCSC
UCSC_LUAD_exprs <- paste0(path,'\\UCSC_LUAD\\LUAD_exprs.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
UCSC_LUAD_exprs=data.frame(t(UCSC_LUAD_exprs))
colnames(UCSC_LUAD_exprs) <- UCSC_LUAD_exprs[1,]
UCSC_LUAD_exprs <- UCSC_LUAD_exprs[-1,]
UCSC_LUAD_exprs$Sample <- rownames(UCSC_LUAD_exprs)

UCSC_LUAD <- paste0(path,'\\UCSC_LUAD\\LUAD.txt') %>% 
  read.table(sep="\t",header=T,check.names=F)
rownames(UCSC_LUAD) <- UCSC_LUAD[,1]
UCSC_LUAD_os <- UCSC_LUAD[,c('Sample','os_time','os_status')]

UCSC_LUAD_exprs_os <- merge(UCSC_LUAD_os,UCSC_LUAD_exprs,by = 'Sample')
UCSC_LUAD_exprs_os1 <- UCSC_LUAD_exprs_os[UCSC_LUAD_exprs_os$os_time>30,]
UCSC_LUAD_exprs_os2 <- na.omit(UCSC_LUAD_exprs_os1)
colnames(UCSC_LUAD_exprs_os2)[1:3] <- c('ID','OS.time','OS')

write.table(UCSC_LUAD_exprs_os2, file=paste0(path,'\\UCSC_LUAD\\UCSC_LUAD_exprs_os.csv'),
            sep=",",quote=F,row.names=F)
