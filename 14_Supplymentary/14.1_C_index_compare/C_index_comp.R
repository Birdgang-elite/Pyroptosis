##########################################################
##根据选取的数据集，对代码进行增删即可，无特殊参数需要设置
##源代码中包括GSE30219、GSE30219、GSE41271、GSE50081和UCSC
##########################################################

#GSE30219
rm(list=ls())
path <- getwd()
setwd(path)
library(survival)
library(survminer)
library(ggplot2)
library(compareC)
GSE30219 <- paste0(path,'\\GSE30219\\data_GSE30219.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

# 单因素Cox回归分析并计算C-index及其标准误
cox_fit_Risk_score <- coxph(Surv(os_time, os_status) ~ Risk_score, data = GSE30219)
c_index_Risk_score <- cox_fit_Risk_score$concordance['concordance']
c_index_Risk_score_se <- cox_fit_Risk_score$concordance['std']

cox_fit_Age <- coxph(Surv(os_time, os_status) ~ Age, data = GSE30219)
c_index_Age <- cox_fit_Age$concordance['concordance']
c_index_Age_se <- cox_fit_Age$concordance['std']

cox_fit_Gender <- coxph(Surv(os_time, os_status) ~ Gender, data = GSE30219)
c_index_Gender <- cox_fit_Gender$concordance['concordance']
c_index_Gender_se <- cox_fit_Gender$concordance['std']

cox_fit_pT <- coxph(Surv(os_time, os_status) ~ pT, data = GSE30219)
c_index_pT <- cox_fit_pT$concordance['concordance']
c_index_pT_se <- cox_fit_pT$concordance['std']

cox_fit_pN <- coxph(Surv(os_time, os_status) ~ pN, data = GSE30219)
c_index_pN <- cox_fit_pN$concordance['concordance']
c_index_pN_se <- cox_fit_pN$concordance['std']

cox_fit_Stage <- coxph(Surv(os_time, os_status) ~ Stage, data = GSE30219)
c_index_Stage <- cox_fit_Stage$concordance['concordance']
c_index_Stage_se <- cox_fit_Stage$concordance['std']

# 创建一个数据框用于绘图
c_index_data <- data.frame(
  Variable = c("Risk_score", "Age", "Gender", "pT", "pN", "Stage"),
  C_Index = c(c_index_Risk_score, c_index_Age, c_index_Gender, 
              c_index_pT, c_index_pN, c_index_Stage),
  SE = c(c_index_Risk_score_se, c_index_Age_se, c_index_Gender_se, 
         c_index_pT_se, c_index_pN_se, c_index_Stage_se)
)

# 计算误差条的上下限
c_index_data$Lower <- c_index_data$C_Index - c_index_data$SE
c_index_data$Upper <- c_index_data$C_Index + c_index_data$SE

# 手动进行C-index的统计比较
# 计算z值和p值
z_RS_vs_Risk_score <- (c_index_Risk_score - c_index_Risk_score) / sqrt(c_index_Risk_score_se^2 + c_index_Risk_score_se^2)
z_RS_vs_Age <- (c_index_Risk_score - c_index_Age) / sqrt(c_index_Risk_score_se^2 + c_index_Age_se^2)
z_RS_vs_Gender <- (c_index_Risk_score - c_index_Gender) / sqrt(c_index_Risk_score_se^2 + c_index_Gender_se^2)
z_RS_vs_pT <- (c_index_Risk_score - c_index_pT) / sqrt(c_index_Risk_score_se^2 + c_index_pT_se^2)
z_RS_vs_pN <- (c_index_Risk_score - c_index_pN) / sqrt(c_index_Risk_score_se^2 + c_index_pN_se^2)
z_RS_vs_Stage <- (c_index_Risk_score - c_index_Stage) / sqrt(c_index_Risk_score_se^2 + c_index_Stage_se^2)

p_RS_vs_Risk_score <- 2 * pnorm(-abs(z_RS_vs_Risk_score))
p_RS_vs_Age <- 2 * pnorm(-abs(z_RS_vs_Age))
p_RS_vs_Gender <- 2 * pnorm(-abs(z_RS_vs_Gender))
p_RS_vs_pT <- 2 * pnorm(-abs(z_RS_vs_pT))
p_RS_vs_pN <- 2 * pnorm(-abs(z_RS_vs_pN))
p_RS_vs_Stage <- 2 * pnorm(-abs(z_RS_vs_Stage))

# 定义p值到星号的转换函数
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")  # 不显著
  }
}

#将z、p和统计差异*号添加到c_index_data数据框中输出
c_index_data$z <- c(z_RS_vs_Risk_score, z_RS_vs_Age, z_RS_vs_Gender, z_RS_vs_pT, z_RS_vs_pN, z_RS_vs_Stage)
c_index_data$p <- c(p_RS_vs_Risk_score, p_RS_vs_Age, p_RS_vs_Gender, p_RS_vs_pT, p_RS_vs_pN, p_RS_vs_Stage)
c_index_data$stars <- sapply(c_index_data$p, p_to_stars)
c_index_data$stars[1] <- " "
write.table(c_index_data,paste0(path,"\\GSE30219\\c_index_data_GSE30219.txt"),
            sep="\t",row.names=F,quote=F, na="")

# 绘制柱状图并添加误差条
p1 <- ggplot(c_index_data, aes(x = factor(Variable,levels = c('Risk_score', 'Age', 'Gender', 'pT', 'pN', 'Stage')), y = C_Index, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_bw() +
  geom_text(aes(label = stars), vjust = -1.5) +
  labs(title = "C-Index Comparison of Different Variables",
       x = "Variable",
       y = "C-Index (Compare with Risk score)") +
  scale_fill_brewer(palette = "Pastel1")
ggsave(paste0(path,"\\GSE30219\\GSE30219.pdf"), p1, width=6,height=5)


#GSE31210
rm(list=ls())
path <- getwd()
setwd(path)
library(survival)
library(survminer)
library(ggplot2)
library(compareC)
GSE31210 <- paste0(path,'\\GSE31210\\data_GSE31210.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

# 单因素Cox回归分析并计算C-index及其标准误
cox_fit_Risk_score <- coxph(Surv(os_time, os_status) ~ Risk_score, data = GSE31210)
c_index_Risk_score <- cox_fit_Risk_score$concordance['concordance']
c_index_Risk_score_se <- cox_fit_Risk_score$concordance['std']

cox_fit_Age <- coxph(Surv(os_time, os_status) ~ Age, data = GSE31210)
c_index_Age <- cox_fit_Age$concordance['concordance']
c_index_Age_se <- cox_fit_Age$concordance['std']

cox_fit_Gender <- coxph(Surv(os_time, os_status) ~ Gender, data = GSE31210)
c_index_Gender <- cox_fit_Gender$concordance['concordance']
c_index_Gender_se <- cox_fit_Gender$concordance['std']

cox_fit_Stage <- coxph(Surv(os_time, os_status) ~ Stage, data = GSE31210)
c_index_Stage <- cox_fit_Stage$concordance['concordance']
c_index_Stage_se <- cox_fit_Stage$concordance['std']

cox_fit_Smoking <- coxph(Surv(os_time, os_status) ~ Smoking, data = GSE31210)
c_index_Smoking <- cox_fit_Smoking$concordance['concordance']
c_index_Smoking_se <- cox_fit_Smoking$concordance['std']

# 创建一个数据框用于绘图
c_index_data <- data.frame(
  Variable = c("Risk_score", "Age", "Gender", "Stage", "Smoking"),
  C_Index = c(c_index_Risk_score, c_index_Age, c_index_Gender, 
              c_index_Stage,c_index_Smoking),
  SE = c(c_index_Risk_score_se, c_index_Age_se, c_index_Gender_se, 
         c_index_Stage_se,c_index_Smoking_se)
)

# 计算误差条的上下限
c_index_data$Lower <- c_index_data$C_Index - c_index_data$SE
c_index_data$Upper <- c_index_data$C_Index + c_index_data$SE

# 手动进行C-index的统计比较
# 计算z值和p值
z_RS_vs_Risk_score <- (c_index_Risk_score - c_index_Risk_score) / sqrt(c_index_Risk_score_se^2 + c_index_Risk_score_se^2)
z_RS_vs_Age <- (c_index_Risk_score - c_index_Age) / sqrt(c_index_Risk_score_se^2 + c_index_Age_se^2)
z_RS_vs_Gender <- (c_index_Risk_score - c_index_Gender) / sqrt(c_index_Risk_score_se^2 + c_index_Gender_se^2)
z_RS_vs_Stage <- (c_index_Risk_score - c_index_Stage) / sqrt(c_index_Risk_score_se^2 + c_index_Stage_se^2)
z_RS_vs_Smoking <- (c_index_Risk_score - c_index_Smoking) / sqrt(c_index_Risk_score_se^2 + c_index_Smoking_se^2)

p_RS_vs_Risk_score <- 2 * pnorm(-abs(z_RS_vs_Risk_score))
p_RS_vs_Age <- 2 * pnorm(-abs(z_RS_vs_Age))
p_RS_vs_Gender <- 2 * pnorm(-abs(z_RS_vs_Gender))
p_RS_vs_Stage <- 2 * pnorm(-abs(z_RS_vs_Stage))
p_RS_vs_Smoking <- 2 * pnorm(-abs(z_RS_vs_Smoking))

# 定义p值到星号的转换函数
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")  # 不显著
  }
}

#将z、p和统计差异*号添加到c_index_data数据框中输出
c_index_data$z <- c(z_RS_vs_Risk_score, z_RS_vs_Age, z_RS_vs_Gender, z_RS_vs_Stage, z_RS_vs_Smoking)
c_index_data$p <- c(p_RS_vs_Risk_score, p_RS_vs_Age, p_RS_vs_Gender, p_RS_vs_Stage, p_RS_vs_Smoking)
c_index_data$stars <- sapply(c_index_data$p, p_to_stars)
c_index_data$stars[1] <- " "
write.table(c_index_data,paste0(path,"\\GSE31210\\c_index_data_GSE31210.txt"),
            sep="\t",row.names=F,quote=F, na="")

# 绘制柱状图并添加误差条
p1 <- ggplot(c_index_data, aes(x = factor(Variable,levels = c('Risk_score', 'Age', 'Gender', 'Stage', 'Smoking')), y = C_Index, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_bw() +
  geom_text(aes(label = stars), vjust = -2.5) +
  labs(title = "C-Index Comparison of Different Variables",
       x = "Variable",
       y = "C-Index (Compare with Risk score)") +
  scale_fill_brewer(palette = "Pastel1")
ggsave(paste0(path,"\\GSE31210\\GSE31210.pdf"), p1, width=6,height=5)

#GSE41271
rm(list=ls())
path <- getwd()
setwd(path)
library(survival)
library(survminer)
library(ggplot2)
library(compareC)
GSE41271 <- paste0(path,'\\GSE41271\\data_GSE41271.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

# 单因素Cox回归分析并计算C-index及其标准误
cox_fit_Risk_score <- coxph(Surv(os_time, os_status) ~ Risk_score, data = GSE41271)
c_index_Risk_score <- cox_fit_Risk_score$concordance['concordance']
c_index_Risk_score_se <- cox_fit_Risk_score$concordance['std']

cox_fit_Age <- coxph(Surv(os_time, os_status) ~ Age, data = GSE41271)
c_index_Age <- cox_fit_Age$concordance['concordance']
c_index_Age_se <- cox_fit_Age$concordance['std']

cox_fit_Gender <- coxph(Surv(os_time, os_status) ~ Gender, data = GSE41271)
c_index_Gender <- cox_fit_Gender$concordance['concordance']
c_index_Gender_se <- cox_fit_Gender$concordance['std']

cox_fit_Stage <- coxph(Surv(os_time, os_status) ~ Stage, data = GSE41271)
c_index_Stage <- cox_fit_Stage$concordance['concordance']
c_index_Stage_se <- cox_fit_Stage$concordance['std']

cox_fit_Smoking <- coxph(Surv(os_time, os_status) ~ Smoking, data = GSE41271)
c_index_Smoking <- cox_fit_Smoking$concordance['concordance']
c_index_Smoking_se <- cox_fit_Smoking$concordance['std']

# 创建一个数据框用于绘图
c_index_data <- data.frame(
  Variable = c("Risk_score", "Age", "Gender", "Stage", "Smoking"),
  C_Index = c(c_index_Risk_score, c_index_Age, c_index_Gender, 
              c_index_Stage,c_index_Smoking),
  SE = c(c_index_Risk_score_se, c_index_Age_se, c_index_Gender_se, 
         c_index_Stage_se,c_index_Smoking_se)
)

# 计算误差条的上下限
c_index_data$Lower <- c_index_data$C_Index - c_index_data$SE
c_index_data$Upper <- c_index_data$C_Index + c_index_data$SE

# 手动进行C-index的统计比较
# 计算z值和p值
z_RS_vs_Risk_score <- (c_index_Risk_score - c_index_Risk_score) / sqrt(c_index_Risk_score_se^2 + c_index_Risk_score_se^2)
z_RS_vs_Age <- (c_index_Risk_score - c_index_Age) / sqrt(c_index_Risk_score_se^2 + c_index_Age_se^2)
z_RS_vs_Gender <- (c_index_Risk_score - c_index_Gender) / sqrt(c_index_Risk_score_se^2 + c_index_Gender_se^2)
z_RS_vs_Stage <- (c_index_Risk_score - c_index_Stage) / sqrt(c_index_Risk_score_se^2 + c_index_Stage_se^2)
z_RS_vs_Smoking <- (c_index_Risk_score - c_index_Smoking) / sqrt(c_index_Risk_score_se^2 + c_index_Smoking_se^2)

p_RS_vs_Risk_score <- 2 * pnorm(-abs(z_RS_vs_Risk_score))
p_RS_vs_Age <- 2 * pnorm(-abs(z_RS_vs_Age))
p_RS_vs_Gender <- 2 * pnorm(-abs(z_RS_vs_Gender))
p_RS_vs_Stage <- 2 * pnorm(-abs(z_RS_vs_Stage))
p_RS_vs_Smoking <- 2 * pnorm(-abs(z_RS_vs_Smoking))

# 定义p值到星号的转换函数
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")  # 不显著
  }
}

#将z、p和统计差异*号添加到c_index_data数据框中输出
c_index_data$z <- c(z_RS_vs_Risk_score, z_RS_vs_Age, z_RS_vs_Gender, z_RS_vs_Stage, z_RS_vs_Smoking)
c_index_data$p <- c(p_RS_vs_Risk_score, p_RS_vs_Age, p_RS_vs_Gender, p_RS_vs_Stage, p_RS_vs_Smoking)
c_index_data$stars <- sapply(c_index_data$p, p_to_stars)
c_index_data$stars[1] <- " "
write.table(c_index_data,paste0(path,"\\GSE41271\\c_index_data_GSE41271.txt"),
            sep="\t",row.names=F,quote=F, na="")

# 绘制柱状图并添加误差条
p1 <- ggplot(c_index_data, aes(x = factor(Variable,levels = c('Risk_score', 'Age', 'Gender', 'Stage', 'Smoking')), y = C_Index, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_bw() +
  geom_text(aes(label = stars), vjust = -1.8) +
  labs(title = "C-Index Comparison of Different Variables",
       x = "Variable",
       y = "C-Index (Compare with Risk score)") +
  scale_fill_brewer(palette = "Pastel1")
ggsave(paste0(path,"\\GSE41271\\GSE41271.pdf"), p1, width=6,height=5)

#GSE50081
rm(list=ls())
path <- getwd()
setwd(path)
library(survival)
library(survminer)
library(ggplot2)
library(compareC)
GSE50081 <- paste0(path,'\\GSE50081\\data_GSE50081.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

# 单因素Cox回归分析并计算C-index及其标准误
cox_fit_Risk_score <- coxph(Surv(os_time, os_status) ~ Risk_score, data = GSE50081)
c_index_Risk_score <- cox_fit_Risk_score$concordance['concordance']
c_index_Risk_score_se <- cox_fit_Risk_score$concordance['std']

cox_fit_Age <- coxph(Surv(os_time, os_status) ~ Age, data = GSE50081)
c_index_Age <- cox_fit_Age$concordance['concordance']
c_index_Age_se <- cox_fit_Age$concordance['std']

cox_fit_Gender <- coxph(Surv(os_time, os_status) ~ Gender, data = GSE50081)
c_index_Gender <- cox_fit_Gender$concordance['concordance']
c_index_Gender_se <- cox_fit_Gender$concordance['std']

cox_fit_pT <- coxph(Surv(os_time, os_status) ~ pT, data = GSE50081)
c_index_pT <- cox_fit_pT$concordance['concordance']
c_index_pT_se <- cox_fit_pT$concordance['std']

cox_fit_pN <- coxph(Surv(os_time, os_status) ~ pN, data = GSE50081)
c_index_pN <- cox_fit_pN$concordance['concordance']
c_index_pN_se <- cox_fit_pN$concordance['std']

cox_fit_Stage <- coxph(Surv(os_time, os_status) ~ Stage, data = GSE50081)
c_index_Stage <- cox_fit_Stage$concordance['concordance']
c_index_Stage_se <- cox_fit_Stage$concordance['std']

cox_fit_Smoking <- coxph(Surv(os_time, os_status) ~ Smoking, data = GSE50081)
c_index_Smoking <- cox_fit_Smoking$concordance['concordance']
c_index_Smoking_se <- cox_fit_Smoking$concordance['std']


# 创建一个数据框用于绘图
c_index_data <- data.frame(
  Variable = c("Risk_score", "Age", "Gender", "pT", "pN", "Stage", "Smoking"),
  C_Index = c(c_index_Risk_score, c_index_Age, c_index_Gender, 
              c_index_pT, c_index_pN, c_index_Stage, c_index_Smoking),
  SE = c(c_index_Risk_score_se, c_index_Age_se, c_index_Gender_se, 
         c_index_pT_se, c_index_pN_se, c_index_Stage_se, c_index_Smoking_se)
)

# 计算误差条的上下限
c_index_data$Lower <- c_index_data$C_Index - c_index_data$SE
c_index_data$Upper <- c_index_data$C_Index + c_index_data$SE

# 手动进行C-index的统计比较
# 计算z值和p值
z_RS_vs_Risk_score <- (c_index_Risk_score - c_index_Risk_score) / sqrt(c_index_Risk_score_se^2 + c_index_Risk_score_se^2)
z_RS_vs_Age <- (c_index_Risk_score - c_index_Age) / sqrt(c_index_Risk_score_se^2 + c_index_Age_se^2)
z_RS_vs_Gender <- (c_index_Risk_score - c_index_Gender) / sqrt(c_index_Risk_score_se^2 + c_index_Gender_se^2)
z_RS_vs_pT <- (c_index_Risk_score - c_index_pT) / sqrt(c_index_Risk_score_se^2 + c_index_pT_se^2)
z_RS_vs_pN <- (c_index_Risk_score - c_index_pN) / sqrt(c_index_Risk_score_se^2 + c_index_pN_se^2)
z_RS_vs_Stage <- (c_index_Risk_score - c_index_Stage) / sqrt(c_index_Risk_score_se^2 + c_index_Stage_se^2)
z_RS_vs_Smoking <- (c_index_Risk_score - c_index_Smoking) / sqrt(c_index_Risk_score_se^2 + c_index_Smoking_se^2)

p_RS_vs_Risk_score <- 2 * pnorm(-abs(z_RS_vs_Risk_score))
p_RS_vs_Age <- 2 * pnorm(-abs(z_RS_vs_Age))
p_RS_vs_Gender <- 2 * pnorm(-abs(z_RS_vs_Gender))
p_RS_vs_pT <- 2 * pnorm(-abs(z_RS_vs_pT))
p_RS_vs_pN <- 2 * pnorm(-abs(z_RS_vs_pN))
p_RS_vs_Stage <- 2 * pnorm(-abs(z_RS_vs_Stage))
p_RS_vs_Smoking <- 2 * pnorm(-abs(z_RS_vs_Smoking))


# 定义p值到星号的转换函数
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")  # 不显著
  }
}


#将z、p和统计差异*号添加到c_index_data数据框中输出
c_index_data$z <- c(z_RS_vs_Risk_score, z_RS_vs_Age, z_RS_vs_Gender, z_RS_vs_pT, z_RS_vs_pN, z_RS_vs_Stage, z_RS_vs_Smoking)
c_index_data$p <- c(p_RS_vs_Risk_score, p_RS_vs_Age, p_RS_vs_Gender, p_RS_vs_pT, p_RS_vs_pN, p_RS_vs_Stage, p_RS_vs_Smoking)
c_index_data$stars <- sapply(c_index_data$p, p_to_stars)
c_index_data$stars[1] <- " "
write.table(c_index_data,paste0(path,"\\GSE50081\\c_index_data_GSE50081.txt"),
            sep="\t",row.names=F,quote=F, na="")

# 绘制柱状图并添加误差条
p1 <- ggplot(c_index_data, aes(x = factor(Variable,levels = c('Risk_score', 'Age', 'Gender', 'pT', 'pN', 'Stage', 'Smoking')), y = C_Index, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_bw() +
  geom_text(aes(label = stars), vjust = -1.5) +
  labs(title = "C-Index Comparison of Different Variables",
       x = "Variable",
       y = "C-Index (Compare with Risk score)") +
  scale_fill_brewer(palette = "Pastel1")
ggsave(paste0(path,"\\GSE50081\\GSE50081.pdf"), p1, width=6,height=5)

#UCSC
rm(list=ls())
path <- getwd()
setwd(path)
library(survival)
library(survminer)
library(ggplot2)
library(compareC)
UCSC <- paste0(path,'\\UCSC_LUAD\\data_all.txt') %>%
  read.table(sep="\t",header=T,check.names=F)

# 单因素Cox回归分析并计算C-index及其标准误
cox_fit_Risk_score <- coxph(Surv(os_time, os_status) ~ Risk_score, data = UCSC)
c_index_Risk_score <- cox_fit_Risk_score$concordance['concordance']
c_index_Risk_score_se <- cox_fit_Risk_score$concordance['std']

cox_fit_Age <- coxph(Surv(os_time, os_status) ~ Age, data = UCSC)
c_index_Age <- cox_fit_Age$concordance['concordance']
c_index_Age_se <- cox_fit_Age$concordance['std']

cox_fit_Gender <- coxph(Surv(os_time, os_status) ~ Gender, data = UCSC)
c_index_Gender <- cox_fit_Gender$concordance['concordance']
c_index_Gender_se <- cox_fit_Gender$concordance['std']

cox_fit_pT <- coxph(Surv(os_time, os_status) ~ pT, data = UCSC)
c_index_pT <- cox_fit_pT$concordance['concordance']
c_index_pT_se <- cox_fit_pT$concordance['std']

cox_fit_pN <- coxph(Surv(os_time, os_status) ~ pN, data = UCSC)
c_index_pN <- cox_fit_pN$concordance['concordance']
c_index_pN_se <- cox_fit_pN$concordance['std']

cox_fit_Stage <- coxph(Surv(os_time, os_status) ~ Stage, data = UCSC)
c_index_Stage <- cox_fit_Stage$concordance['concordance']
c_index_Stage_se <- cox_fit_Stage$concordance['std']

cox_fit_Smoking <- coxph(Surv(os_time, os_status) ~ Smoking, data = UCSC)
c_index_Smoking <- cox_fit_Smoking$concordance['concordance']
c_index_Smoking_se <- cox_fit_Smoking$concordance['std']


# 创建一个数据框用于绘图
c_index_data <- data.frame(
  Variable = c("Risk_score", "Age", "Gender", "pT", "pN", "Stage", "Smoking"),
  C_Index = c(c_index_Risk_score, c_index_Age, c_index_Gender, 
              c_index_pT, c_index_pN, c_index_Stage, c_index_Smoking),
  SE = c(c_index_Risk_score_se, c_index_Age_se, c_index_Gender_se, 
         c_index_pT_se, c_index_pN_se, c_index_Stage_se, c_index_Smoking_se)
)

# 计算误差条的上下限
c_index_data$Lower <- c_index_data$C_Index - c_index_data$SE
c_index_data$Upper <- c_index_data$C_Index + c_index_data$SE

# 手动进行C-index的统计比较
# 计算z值和p值
z_RS_vs_Risk_score <- (c_index_Risk_score - c_index_Risk_score) / sqrt(c_index_Risk_score_se^2 + c_index_Risk_score_se^2)
z_RS_vs_Age <- (c_index_Risk_score - c_index_Age) / sqrt(c_index_Risk_score_se^2 + c_index_Age_se^2)
z_RS_vs_Gender <- (c_index_Risk_score - c_index_Gender) / sqrt(c_index_Risk_score_se^2 + c_index_Gender_se^2)
z_RS_vs_pT <- (c_index_Risk_score - c_index_pT) / sqrt(c_index_Risk_score_se^2 + c_index_pT_se^2)
z_RS_vs_pN <- (c_index_Risk_score - c_index_pN) / sqrt(c_index_Risk_score_se^2 + c_index_pN_se^2)
z_RS_vs_Stage <- (c_index_Risk_score - c_index_Stage) / sqrt(c_index_Risk_score_se^2 + c_index_Stage_se^2)
z_RS_vs_Smoking <- (c_index_Risk_score - c_index_Smoking) / sqrt(c_index_Risk_score_se^2 + c_index_Smoking_se^2)

p_RS_vs_Risk_score <- 2 * pnorm(-abs(z_RS_vs_Risk_score))
p_RS_vs_Age <- 2 * pnorm(-abs(z_RS_vs_Age))
p_RS_vs_Gender <- 2 * pnorm(-abs(z_RS_vs_Gender))
p_RS_vs_pT <- 2 * pnorm(-abs(z_RS_vs_pT))
p_RS_vs_pN <- 2 * pnorm(-abs(z_RS_vs_pN))
p_RS_vs_Stage <- 2 * pnorm(-abs(z_RS_vs_Stage))
p_RS_vs_Smoking <- 2 * pnorm(-abs(z_RS_vs_Smoking))


# 定义p值到星号的转换函数
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")  # 不显著
  }
}


#将z、p和统计差异*号添加到c_index_data数据框中输出
c_index_data$z <- c(z_RS_vs_Risk_score, z_RS_vs_Age, z_RS_vs_Gender, z_RS_vs_pT, z_RS_vs_pN, z_RS_vs_Stage, z_RS_vs_Smoking)
c_index_data$p <- c(p_RS_vs_Risk_score, p_RS_vs_Age, p_RS_vs_Gender, p_RS_vs_pT, p_RS_vs_pN, p_RS_vs_Stage, p_RS_vs_Smoking)
c_index_data$stars <- sapply(c_index_data$p, p_to_stars)
c_index_data$stars[1] <- " "
write.table(c_index_data,paste0(path,"\\UCSC_LUAD\\c_index_data_UCSC.txt"),
            sep="\t",row.names=F,quote=F, na="")

# 绘制柱状图并添加误差条
p1 <- ggplot(c_index_data, aes(x = factor(Variable,levels = c('Risk_score', 'Age', 'Gender', 'pT', 'pN', 'Stage', 'Smoking')), y = C_Index, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_bw() +
  geom_text(aes(label = stars), vjust = -1.5) +
  labs(title = "C-Index Comparison of Different Variables",
       x = "Variable",
       y = "C-Index (Compare with Risk score)") +
  scale_fill_brewer(palette = "Pastel1")
ggsave(paste0(path,"\\UCSC_LUAD\\UCSC.pdf"), p1, width=6,height=5)