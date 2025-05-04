###########################################
##基于RS、T、N构建列线图、Calibration和DCA。
###########################################

rm(list=ls())
path <- getwd()
setwd(path) 
library(survival)
library(rms)

#Data_preparing
data <- read.table("data_all.txt",sep="\t",header=T,check.names=F) 
cli_exp <- data
ddist <- datadist(cli_exp)
options(datadist='ddist')


#Nomogram 
cox <- cph(Surv(os_time,os_status) ~ pT+pN+Risk_score,surv=T,x=T, y=T,data=cli_exp)
surv <- Survival(cox)
sur_1_year<-function(x)surv(365,lp=x)#1年生存
sur_3_year<-function(x)surv(1095,lp=x)#3年生存
sur_5_year<-function(x)surv(1825,lp=x)#5年生存

nom_sur<- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,
                   funlabel=c('1-Year Survival','3-Year Survival','5-Year Survival'),
                   maxscale=100,
                   fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3'))
pdf("Nomogram.pdf",6,4)
plot(nom_sur)
dev.off()

#Calibration curve
rm(list=ls())
path <- getwd()
setwd(path) 
cli_exp <-read.table("data_all.txt",sep="\t",header=T,check.names=F)
library("survival")
library("survminer")
library(rms)

ddist <- datadist(cli_exp)
options(datadist='ddist')
cli_exp[,"os_time"]=cli_exp[,"os_time"]/365  #时间要转换为年
units(cli_exp$os_time) <- "Year"

#####1 year calibration curve#####
cox1 <- cph(Surv(os_time,os_status) ~ pT+pN+Risk_score,surv=T,x=T, y=T,time.inc = 1,data=cli_exp)
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u=1, m= 100, B=1000)
pdf("calibrate_os1.pdf")
plot(cal1,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Risk_score-Predicted Probability of 1-Year Survival",ylab="Actual 1-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()

#####3 year calibration curve#####
cox3 <- cph(Surv(os_time,os_status) ~ pT+pN+Risk_score,surv=T,x=T, y=T,time.inc = 3,data=cli_exp)
cal3 <- calibrate(cox3, cmethod="KM", method="boot", u=3, m= 100, B=1000)
pdf("calibrate_os3.pdf")
plot(cal3,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Risk_score-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()

#####5 year calibration curve#####
cox5 <- cph(Surv(os_time,os_status) ~ pT+pN+Risk_score,surv=T,x=T, y=T,time.inc = 5,data=cli_exp)
cal5 <- calibrate(cox5, cmethod="KM", method="boot", u=5, m= 100, B=1000)
pdf("calibrate_os5.pdf")
plot(cal5,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Risk_score-Predicted Probability of 5-Year Survival",ylab="Actual 5-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()


#1 year DCA curve
rm(list=ls())
path <- getwd()
setwd(path) 
library('dcurves')
library('survival')
library('ggplot2')
library('dplyr')
data.set <- read.table("data_all.txt",sep="\t",header=T,check.names=F)

data.set = data.set[,c('Sample','os_status','os_time','pT','pN','Risk_score')]
data.set = na.omit(data.set)

Srv = Surv(data.set$os_time, data.set$os_status)


coxmod2 = coxph(Srv ~ pT, data=data.set)
coxmod3 = coxph(Srv ~ pN, data=data.set)
coxmod4 = coxph(Srv ~ Risk_score, data=data.set)
coxmod5 = coxph(Srv ~ pT+pN+Risk_score, data=data.set)



data.set$pT1 = c(1- (summary(survfit(coxmod2,newdata=data.set), times=365)$surv))
data.set$pN1 = c(1- (summary(survfit(coxmod3,newdata=data.set), times=365)$surv))
data.set$Risk_score1 = c(1- (summary(survfit(coxmod4,newdata=data.set), times=365)$surv))
data.set$Nom = c(1- (summary(survfit(coxmod5,newdata=data.set), times=365)$surv))

pdf("DCA1.pdf",5,5)
dca(Srv~pT1+pN1+Risk_score1+Nom,data = data.set,time = 365) %>%plot(show_ggplot_code = T)
dev.off()



#3 year DCA curve
rm(list=ls())
path <- getwd()
setwd(path) 
library('dcurves')
library('survival')
library('ggplot2')
library('dplyr')
data.set <- read.table("data_all.txt",sep="\t",header=T,check.names=F)

data.set = data.set[,c('Sample','os_status','os_time','pT','pN','Risk_score')]
data.set = na.omit(data.set)

Srv = Surv(data.set$os_time, data.set$os_status)


coxmod2 = coxph(Srv ~ pT, data=data.set)
coxmod3 = coxph(Srv ~ pN, data=data.set)
coxmod4 = coxph(Srv ~ Risk_score, data=data.set)
coxmod5 = coxph(Srv ~ pT+pN+Risk_score, data=data.set)



data.set$pT3 = c(1- (summary(survfit(coxmod2,newdata=data.set), times=1095)$surv))
data.set$pN3 = c(1- (summary(survfit(coxmod3,newdata=data.set), times=1095)$surv))
data.set$Risk_score3 = c(1- (summary(survfit(coxmod4,newdata=data.set), times=1095)$surv))
data.set$Nom = c(1- (summary(survfit(coxmod5,newdata=data.set), times=1095)$surv))

pdf("DCA3.pdf",5,5)
dca(Srv~pT3+pN3+Risk_score3+Nom,data = data.set,time = 1095) %>%plot(show_ggplot_code = T)
dev.off()


#5 year DCA curve
rm(list=ls())
path <- getwd()
setwd(path) 
library('dcurves')
library('survival')
library('ggplot2')
library('dplyr')
data.set <- read.table("data_all.txt",sep="\t",header=T,check.names=F)

data.set = data.set[,c('Sample','os_status','os_time','pT','pN','Risk_score')]
data.set = na.omit(data.set)

Srv = Surv(data.set$os_time, data.set$os_status)


coxmod2 = coxph(Srv ~ pT, data=data.set)
coxmod3 = coxph(Srv ~ pN, data=data.set)
coxmod4 = coxph(Srv ~ Risk_score, data=data.set)
coxmod5 = coxph(Srv ~ pT+pN+Risk_score, data=data.set)



data.set$pT5 = c(1- (summary(survfit(coxmod2,newdata=data.set), times=1825)$surv))
data.set$pN5 = c(1- (summary(survfit(coxmod3,newdata=data.set), times=1825)$surv))
data.set$Risk_score5 = c(1- (summary(survfit(coxmod4,newdata=data.set), times=1825)$surv))
data.set$Nom = c(1- (summary(survfit(coxmod5,newdata=data.set), times=1825)$surv))

pdf("DCA5.pdf",5,5)
dca(Srv~pT5+pN5+Risk_score5+Nom,data = data.set,time = 1825) %>%plot(show_ggplot_code = T)
dev.off()

