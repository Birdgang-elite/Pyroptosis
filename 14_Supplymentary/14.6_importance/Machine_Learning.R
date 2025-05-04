rm(list=ls())
path <- getwd()
setwd(path)

library(Mime1)
library(fastAdaboost)
library(CoxBoost)
library(tidyverse)

#读取表达谱和预后数据
load('res.RData')


message("---1-7.RSF + Ridge ---")
result <- data.frame()
ml.res = list()
riskscore = list()
seed <- 5201314
rf_nodesize <- 5
train_data <- list_train_vali_Data$TCGA
est_dd <- as.data.frame(train_data)[, c('ID','OS.time','OS',res$Sig.genes)]
rownames(est_dd) <- est_dd[,1]
est_dd <- est_dd[,-1]
colnames(est_dd) <- gsub("-","_",colnames(est_dd))


library(randomForestSRC)
library(glmnet)

fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
             importance = T, proximity = T, forest = T, seed = seed)

#绘制贡献度图
############################################
library(ggplot2)
library(dplyr)
library(tibble)
library(viridis)

importance_gene <- data.frame(fit$importance) %>% 
  rownames_to_column("gene") %>% 
  arrange(- fit.importance)
importance_gene

ggplot(data=importance_gene, aes(x = reorder(gene,  fit.importance), 
                                 y=fit.importance,fill=gene)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  scale_fill_viridis(discrete = TRUE, option = "D") +
  #scale_fill_brewer(palette = "Set3") +
  coord_flip()
ggsave("model_genes.pdf",width = 9,height = 7)
#################################################



rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

returnIDtoRS = function(rs.table.list, rawtableID) {
  for (i in names(rs.table.list)) {
    rs.table.list[[i]]$ID = rawtableID[[i]]$ID
    rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select("ID", 
                                                              everything())
  }
  return(rs.table.list)
}

if (length(rid) > 1) {
  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
  val_dd_list2 <- lapply(list_train_vali_Data, 
                         function(x) {
                           x[, c("OS.time", "OS", rid)]
                         })
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", 
                  alpha = 0, type.measure = "class")
  
  
  # 获取最优的lambda值
  best_lambda <- fit$lambda.min
  
  # 使用最优的lambda值拟合最终模型
  final_fit <- glmnet(x1, x2, family = "cox", alpha = 0, lambda = best_lambda)
  
  # 输出最终的特征和系数
  coefficients <- coef(final_fit)
  non_zero_coef <- coefficients[coefficients != 0]
  print("最终的特征和系数：")
  print(non_zero_coef)
  
  
  
  
  
  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, 
                                            type = "response", newx = as.matrix(x[, 
                                                                                  -c(1, 2)]), s = fit$lambda.min)))
  })
  cc <- data.frame(Cindex = sapply(rs, function(x) {
    as.numeric(summary(coxph(Surv(OS.time, OS) ~ 
                               RS, x))$concordance[1])
  })) %>% rownames_to_column("ID")
  cc$Model <- paste0("RSF + ", "Ridge")
  result <- rbind(result, cc)
  ml.res[[paste0("RSF + ", "Ridge")]] = fit
  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
  riskscore[[paste0("RSF + ", "Ridge")]] = rs
}  else{
  warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
}  

