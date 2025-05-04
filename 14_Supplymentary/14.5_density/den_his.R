rm(list=ls())
path <- getwd()
setwd(path)
library('ggplot2')
library('patchwork')
library(tidyverse)

#GSE3141
GSE3141 <- paste0(path,'\\GSE3141\\GSE3141_exprs_RS_test.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

GSE3141_median <- median(GSE3141$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\GSE3141\\his_GSE3141.pdf'),width = 5, height = 5)
hist(GSE3141$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = GSE3141_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(GSE3141_median, 0, 
     labels = paste("Median: ", round(GSE3141_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\GSE3141\\den1_GSE3141.pdf'),width = 5, height = 5)
plot(density(GSE3141$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = GSE3141_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(GSE3141_median, max(density(GSE3141$Risk_score)$y), 
     labels = paste("Median: ", round(GSE3141_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\GSE3141\\den2_GSE3141.pdf'),width = 5, height = 5)
ggplot(GSE3141, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()

#GSE30219
GSE30219 <- paste0(path,'\\GSE30219\\GSE30219_exprs_RS_test.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

GSE30219_median <- median(GSE30219$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\GSE30219\\his_GSE30219.pdf'),width = 5, height = 5)
hist(GSE30219$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = GSE30219_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(GSE30219_median, 0, 
     labels = paste("Median: ", round(GSE30219_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\GSE30219\\den1_GSE30219.pdf'),width = 5, height = 5)
plot(density(GSE30219$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = GSE30219_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(GSE30219_median, max(density(GSE30219$Risk_score)$y), 
     labels = paste("Median: ", round(GSE30219_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\GSE30219\\den2_GSE30219.pdf'),width = 5, height = 5)
ggplot(GSE30219, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()


#GSE31210
GSE31210 <- paste0(path,'\\GSE31210\\GSE31210_exprs_RS_test.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

GSE31210_median <- median(GSE31210$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\GSE31210\\his_GSE31210.pdf'),width = 5, height = 5)
hist(GSE31210$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = GSE31210_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(GSE31210_median, 0, 
     labels = paste("Median: ", round(GSE31210_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\GSE31210\\den1_GSE31210.pdf'),width = 5, height = 5)
plot(density(GSE31210$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = GSE31210_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(GSE31210_median, max(density(GSE31210$Risk_score)$y), 
     labels = paste("Median: ", round(GSE31210_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\GSE31210\\den2_GSE31210.pdf'),width = 5, height = 5)
ggplot(GSE31210, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()


#GSE41271
GSE41271 <- paste0(path,'\\GSE41271\\GSE41271_exprs_RS_test.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

GSE41271_median <- median(GSE41271$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\GSE41271\\his_GSE41271.pdf'),width = 5, height = 5)
hist(GSE41271$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = GSE41271_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(GSE41271_median, 0, 
     labels = paste("Median: ", round(GSE41271_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\GSE41271\\den1_GSE41271.pdf'),width = 5, height = 5)
plot(density(GSE41271$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = GSE41271_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(GSE41271_median, max(density(GSE41271$Risk_score)$y), 
     labels = paste("Median: ", round(GSE41271_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\GSE41271\\den2_GSE41271.pdf'),width = 5, height = 5)
ggplot(GSE41271, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()


#GSE50081
GSE50081 <- paste0(path,'\\GSE50081\\GSE50081_exprs_RS_test.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

GSE50081_median <- median(GSE50081$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\GSE50081\\his_GSE50081.pdf'),width = 5, height = 5)
hist(GSE50081$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = GSE50081_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(GSE50081_median, 0, 
     labels = paste("Median: ", round(GSE50081_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\GSE50081\\den1_GSE50081.pdf'),width = 5, height = 5)
plot(density(GSE50081$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = GSE50081_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(GSE50081_median, max(density(GSE50081$Risk_score)$y), 
     labels = paste("Median: ", round(GSE50081_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\GSE50081\\den2_GSE50081.pdf'),width = 5, height = 5)
ggplot(GSE50081, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()


#ucsc_all
TCGA <- paste0(path,'\\UCSC_LUAD\\UCSC_LUAD_exprs_RS_all.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

TCGA_median <- median(TCGA$Risk_score)

# 绘制基因表达值的直方图
pdf(paste0(path,'\\UCSC_LUAD\\his_TCGA.pdf'),width = 5, height = 5)
hist(TCGA$Risk_score, main = "",
     xlab = "Risk score", 
     col = "lightblue", 
     border = "black") 
# 添加中位数垂直线
abline(v = TCGA_median, col = "red", lwd = 2) 
# 添加中位数文本注释
text(TCGA_median, 0, 
     labels = paste("Median: ", round(TCGA_median, 2)), 
     pos = 4, col = "red")
dev.off()

# 绘制基因表达值的密度图
pdf(paste0(path,'\\UCSC_LUAD\\den1_TCGA.pdf'),width = 5, height = 5)
plot(density(TCGA$Risk_score), main = "",
     xlab = "Risk score", 
     col = "red", 
     lwd = 2)
# 在密度图中添加中位数垂直线
abline(v = TCGA_median, col = "blue", lwd = 2, lty = 2)
# 在密度图中添加中位数文本注释
text(TCGA_median, max(density(TCGA$Risk_score)$y), 
     labels = paste("Median: ", round(TCGA_median, 2)), 
     pos = 4, col = "blue")
dev.off()

# 按分组绘制密度图
pdf(paste0(path,'\\UCSC_LUAD\\den2_TCGA.pdf'),width = 5, height = 5)
ggplot(TCGA, aes(x = Risk_score, fill = factor(Risk_score_LH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Risk score",
       y = "density",
       fill = "Group") +
  theme_minimal()
dev.off()