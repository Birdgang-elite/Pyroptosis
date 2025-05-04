rm(list=ls())
path <- getwd()
setwd(path)

library(Mime1)
library(fastAdaboost)
library(CoxBoost)

#读取表达谱和预后数据
TCGA <- read.csv('UCSC_LUAD_exprs_os.csv')
GSE3141 <- read.csv('GSE3141_LUAD_exprs_os.csv')
GSE30219 <- read.csv('GSE30219_LUAD_exprs_os.csv')
GSE31210 <- read.csv('GSE31210_LUAD_exprs_os.csv')
GSE41271 <- read.csv('GSE41271_LUAD_exprs_os.csv')
GSE50081 <- read.csv('GSE50081_LUAD_exprs_os.csv')


list_train_vali_Data <- list(TCGA,GSE3141,GSE30219, GSE31210,GSE41271,GSE50081)
names(list_train_vali_Data) <- c('TCGA','GSE3141','GSE30219','GSE31210','GSE41271','GSE50081')

#读取DEGOIs
genelist <- read.csv('genelist.csv')
genelist <- as.character(genelist$gene)

#构建101种机器学习模型组合
#该包大大降低了学习成本，可以通过ML.Dev.Prog.Sig函数直接构建
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all',nodesize =5,seed = 5201314 )

#计算所有模型在训练集、验证集的C-index,并绘图
pdf(file = 'Cindex_all.pdf',width = 8,height = 12)
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)
dev.off()

#选取最佳的Model
Best_model <- "RSF + Ridge"    ####看前面输出的pdf文件，选取最佳的model（根据Mean_cindex_in_validation的最大值来确定best model）

#计算特定模型在训练集、验证集的C-index,并绘图
pdf(file = 'Cindex_best_model.pdf',width = 8,height = 8)
cindex_dis_select(res,
                  model=Best_model,  
                  order= names(list_train_vali_Data))
dev.off()

#导出
write.table(res[['Cindex.res']],file = 'Cindex.res.txt', sep = '\t', quote = F, col.names = T, row.names = F)
write.table(res[['Sig.genes']],file = 'unicox.txt', sep = '\t', quote = F, col.names = F, row.names = F)

RS <- data.frame()
for(i in seq_along(res[['riskscore']])){
  for(l in c(1:length(list_train_vali_Data))){
    rt = cbind(res[['riskscore']][[i]][[l]][,c(1,4)],
               rep(names(res[['riskscore']][i]), length(rownames(res[['riskscore']][[i]][[l]][,c(1,4)]))))
    colnames(rt)[3] = "model"
    RS = rbind(RS,rt)
  }
}
write.table(RS,file = 'riskscore.txt', sep = '\t', quote = F, col.names = T, row.names = F)

#查看特定模型下各个数据集的KM曲线情况
survplot <- vector("list",6)  #假设有n个数据集，就改为n
for (i in c(1:6)) {     #假设有n个数据集，就改为c(1:n)
  print(survplot[[i]]<-rs_sur(res, model_name = Best_model,
                              dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9))
  )
}
pdf(file = 'KM_best_model.pdf',width = 15,height = 10)
aplot::plot_list(gglist=survplot,ncol=3)  #根据数据集个数，自行设定ncol数值
dev.off()

#计算所有模型的1年、3年、5年的AUC值
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
#可视化所有模型的1年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_1y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.1y,  #如果需要展示3/5年的数据，此处需要修改
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)     #如果需要展示3/5年的数据，此处也需要修改
dev.off()

#可视化所有模型的3年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_3y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.3y,  #如果需要展示3/5年的数据，此处需要修改
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=3)     #如果需要展示3/5年的数据，此处也需要修改
dev.off()

#可视化所有模型的5年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_5y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.5y,  
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=5)     
dev.off()

#可视化特定模型的5年AUC结果
#这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH,所以图不好看，可以自己重新绘图，这里就不输出了
roc_vis(all.auc.5y,   
        model_name = Best_model,  
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=5)   

#可视化特定模型的1年、3年、5年的AUC结果
pdf(file = 'AUC_135y_best_model.pdf',width = 12,height = 4)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name=Best_model,  
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
dev.off()

#最终模型单因素Cox回归的Meta分析
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res, optimal.model = Best_model,   
                                   type = 'categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
pdf(file = ('meta_unicox_vis.pdf'),width = 10, height = 4,onefile = FALSE)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
dev.off()










#计算所有模型在训练集、验证集的AUC,并绘图、导出数据
##计算所有模型的1年、3年、5年的AUC值
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["TCGA"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
##可视化所有模型的1年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_1y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.1y,  #如果需要展示3/5年的数据，此处需要修改
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)     #如果需要展示3/5年的数据，此处也需要修改
dev.off()

#可视化所有模型的3年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_3y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.3y,  #如果需要展示3/5年的数据，此处需要修改
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=3)     #如果需要展示3/5年的数据，此处也需要修改
dev.off()

#可视化所有模型的5年AUC结果（这里使用的是Riskscore绘制的ROC曲线，不是Risk_score_LH）
pdf(file = 'AUC_all_5y.pdf',width = 8,height = 12)
auc_dis_all(all.auc.5y,  
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=5)     
dev.off()

#输出上述数据
##长数据格式
Cindex_all_long <- res[['Cindex.res']]
write.table(Cindex_all_long,file = 'Cindex_all_long.txt', sep = '\t', quote = F, col.names = T, row.names = F)
##宽数据格式
Cindex_all_wide <- spread(Cindex_all_long, key = 'ID', value = 'Cindex')
write.table(Cindex_all_wide,file = 'Cindex_all_wide.txt', sep = '\t', quote = F, col.names = T, row.names = F)

