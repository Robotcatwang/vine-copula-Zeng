## ----setup, echo=F-------------------------------------------------------

################# 第 3 章 R 程序代码  ####################


knitr::opts_knit$set(root.dir = getwd())
knitr::opts_chunk$set(echo = FALSE, results = 'hide')
knitr::opts_chunk$set(warning = FALSE, message=FALSE)


## ----prepare-------------------------------------------------------------
rm(list=ls())
options(digits=4)
options(scipen=100)
graphics.off()
Sys.setlocale("LC_ALL", "Chinese")


## ----import-data---------------------------------------------------------
## 准备
rm(list=ls())
options(digits=4)
options(scipen=100)
graphics.off()
## path <- "F:/github_repo/gas-gold-Wei/" 
## setwd(path)
#Sys.setlocale('LC_ALL','C') 
Sys.setlocale("LC_ALL", "Chinese")


## 加载包
library(zoo)
library(forecast)
library(tseries)
library(xts)
library(ggplot2)
library(FinTS)
library(moments)
library(devtools)
library(FinTS)
library(GAS)
library(fBasics)
library(knitr)
library(timeDate)
library(timeSeries)

## 载入数据
dat=read.csv("./data/jiesuan.csv",header=T)
shibor<-dat[,2]
shibor<-shibor[1:2484]
shibor.date<-as.Date(dat[,1])
shibor.date<-shibor.date[1:2484]

#计算对数收益率
shibor.rt<-diff(log(shibor))
shibor.rt<-shibor.rt[1:2483]
shibor.rt = na.omit(shibor.rt)


## ----stat1, results='markup'---------------------------------------------
shibor.rt1<-basicStats(shibor.rt)
shibor.rt2<- cbind("基本统计量",as.data.frame(t(shibor.rt1[c(7,14,15,16),])))
colnames(shibor.rt2) <- c("黄金期货对数收益率","均值","标准差","偏态系数","峰态系数")

kable(shibor.rt2,row.names =F,align = "c", caption="黄金期货对数收益率描述性分析结果",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----fig-ts1,eval=T,fig.cap = "我国黄金期货主力合约对数收益率时序图", dev='png'----
par(mfrow=c(2,1))
par(mar=c(2,4,1,2))
plot(shibor.date,shibor,type="l",xlab="日期",ylab="每日结算价")
plot(shibor.date[-1],shibor.rt,type="l",xlab="日期",ylab="黄金期货对数收益率",cex.main=0.95,las=1)
#显然具有波动应聚集的特征


## ------------------------------------------------------------------------
#平稳性检验
pp.test(shibor.rt)  #平稳


## ----hist, fig.cap="我国黄金期货主力合约收益率频率分布直方图和正态Q-Q图", dev='png'----
par(mfrow=c(1,2))
#绘制直方图
hist(shibor.rt,main=NULL,breaks = 50)

#正态性检验
qqnorm(shibor.rt,main =NULL)
qqline(shibor.rt)                #QQ图检验非正态


## ------------------------------------------------------------------------
#正态性检验
jarque.bera.test(shibor.rt)    #JB检验非正态


## ------------------------------------------------------------------------
acf(shibor.rt,plot = F,lag.max = 15)
pacf(shibor.rt,plot = F,lag.max = 15)
Box1=Box.test(shibor.rt,lag=10,type = "Ljung")
library(forecast)
library(stats)
fit1=Arima(shibor.rt,order = c(1,0,1))
fit1.res=c("fit1",AIC(fit1),BIC(fit1))
fit2=arima(shibor.rt,order = c(1,0,2))
fit2.res=c("fit2",AIC(fit2),BIC(fit2))
fit3=arima(shibor.rt,order = c(2,0,1))
fit3.res=c("fit3",AIC(fit3),BIC(fit3))
fit4=arima(shibor.rt,order = c(2,0,2))
fit4.res=c("fit4",AIC(fit4),BIC(fit4))
res1=data.frame(fit1.res,fit2.res,fit3.res,fit4.res)
res1
auto.arima(shibor.rt)
fit111=arima(shibor.rt,order = c(1,0,1))
resi<-resid(fit111,standardize=F)
res.ts<-ts(resi,frequency = 250)
Box2=Box.test(resi,lag=10,type = "Ljung")


## ----acf-pacf, fig.cap="我国黄金期货主力合约收益率自相关图和偏自相关图",fig.height=8,dev='png'----
#自相关性检验
par(mfrow=c(2,1)) 
acf(shibor.rt,main="",xlab="滞后期",ylab="ACF")#画自相关图
title(main = "(a)the ACF of Return",cex.main=0.95)
pacf(shibor.rt,main="",xlab="滞后期",ylab="PACF",las=1)#画偏自相关图
title(main="(b)the PACF of Return",cex.main=0.95)


## ------------------------------------------------------------------------
fit=stats::arima(shibor.rt,order = c(1,0,1))
shibor.rt.r<-stats::residuals(fit)


## ----residual, fig.cap="我国黄金期货主力合约收益率残差序列图", dev='png'----
plot(shibor.date[-1],shibor.rt.r,type="l",
     xlab="日期",ylab="残差",cex.main=0.95,las=1)


## ------------------------------------------------------------------------
#ARCH效应检验
ArchTest(shibor.rt.r,lag=12)  #存在ARCH效应


## ----arma-fit, fig.cap="残差平方相关图",fig.height=8,dev="png"-----------
#拟合GARCH
#残差平方的自相关性分析
par(mfrow=c(2,1))
rt.square<-shibor.rt.r^2
acf(rt.square,main="",xlab="lag(c)",ylab="ACF",las=1)#画自相关图
title(main = "(c)the ACF of resi Square",cex.main=0.95)
pacf(rt.square,main="",xlab="Lag(d)",ylab="PACF",las=1)#画偏自相关图
title(main = "(d)the PACF of resi Square",cex.main=0.95)
Box.test(rt.square,lag = 10,type = "Ljung")


## ------------------------------------------------------------------------
gasspec=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=FALSE,scale=FALSE,shape=FALSE))
gasspec0=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=FALSE,scale=TRUE,shape=FALSE))
gasspec1=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=FALSE,scale=TRUE,shape=TRUE))
gasspec2=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=TRUE,scale=TRUE,shape=TRUE))

InSampleData987 = shibor.rt.r[1:1733]
gasfit=UniGASFit(gasspec,data=InSampleData987,fn.optimizer = fn.optim, Compute.SE = TRUE)
gasfit0=UniGASFit(gasspec0,data=InSampleData987,fn.optimizer = fn.optim, Compute.SE = TRUE)
gasfit1=UniGASFit(gasspec1,data=InSampleData987,fn.optimizer = fn.optim, Compute.SE = TRUE)
gasfit2=UniGASFit(gasspec2,data=InSampleData987,fn.optimizer = fn.optim, Compute.SE = TRUE)

L1RT1<-2*(gasfit0@GASDyn$dLLK-gasfit@GASDyn$dLLK)
P1_value1<-1-pchisq(L1RT1,1) 

L1RT2<-2*(gasfit1@GASDyn$dLLK-gasfit0@GASDyn$dLLK)
P1_value2<-1-pchisq(L1RT2,1) 

L1RT3<-2*(gasfit2@GASDyn$dLLK-gasfit1@GASDyn$dLLK)
P1_value3<-1-pchisq(L1RT3,1) 


## ----test-LRT-tab,eval=T,echo=F------------------------------------------
test_LRT=data.frame(a=c(L1RT1,L1RT2,L1RT3),
                    b=c(P1_value1,P1_value2,P1_value3),
                    c=c(1733,1733,1733))
test_LRT<- cbind(c("u1 vs u0","u2 vs u1","u3 vs u2"), test_LRT)
colnames(test_LRT)=c("假设","LRT","P值","样本量")
rownames(test_LRT)=c("u1 vs u0","u2 vs u1","u3 vs u2")
write.csv(test_LRT, "testLRT.csv", row.names = F)


## ----test-LRT, eval=T,results='markup', cache=F--------------------------
tablrt <- read.csv('./testLRT.csv')
knitr::kable(tablrt, row.names =F, align = "l", caption="似然比检验对GAS（1，1)模型的测试结果",
      longtable = TRUE, booktabs = TRUE, linesep  = "", escape = F)


## ----fig8,echo=FALSE,fig.cap="模型估计结果比较表",fig.height=20,fig.width=40,cache=F,dev="png",results='markup'----
knitr::include_graphics("./SELECT1.png")


## ------------------------------------------------------------------------
InSampleData = shibor.rt.r[1:1733]
OutSampleData = shibor.rt.r[1734:2483]

herspec=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=TRUE,scale=TRUE,shape=TRUE))

herfit=UniGASFit(herspec,data=InSampleData,fn.optimizer = fn.optim, Compute.SE = TRUE)


## ----best-model-tab,eval=T,echo=F----------------------------------------
best_model=data.frame(a=c(herfit@Estimates[["lParList"]][["vKappa"]][1],herfit@Estimates[["lParList"]][["vKappa"]][2],herfit@Estimates[["lParList"]][["vKappa"]][3],herfit@Estimates[["lParList"]][["mA"]][1,1],herfit@Estimates[["lParList"]][["mA"]][2,2],herfit@Estimates[["lParList"]][["mA"]][3,3],herfit@Estimates[["lParList"]][["mB"]][1,1],herfit@Estimates[["lParList"]][["mB"]][2,2],herfit@Estimates[["lParList"]][["mB"]][3,3]),
                    b=c(1.10E-04,3.87E-02,1.12E-01,1.19E-09,2.82E-02,3.73E-03,2.66E-05,4.01E-03,6.30E-03),
                    c=c(-0.5593643,-4.867983,-7.7252002,837.5653786,6.7397406,254.4799865,18762.6533,244.1425435,100.4244369),
                    d=c(2.88E-01,5.64E-07,5.55E-15,0.00E+00,7.93E-12,0.00E+00,0.00E+00,0.00E+00,0.00E+00))
best_model <- cbind(c("k1","k2","k3","a1","a2","a3","b1","b2","b3"), best_model)
colnames(best_model)=c("参数","估计值","标准误","t值","P值")
rownames(best_model)=c("k1","k2","k3","a1","a2","a3","b1","b2","b3")
write.csv(best_model, "bestmodel.csv", row.names = F)


## ----best-model, eval=T,results='markup', cache=F------------------------
bestmodel <- read.csv('./bestmodel.csv')
knitr::kable(bestmodel, row.names =F, align = "l", caption="最佳模型估计结果",
      longtable = TRUE, booktabs = TRUE, linesep  = "", escape = F,digits=6)


## ----fig10,echo=FALSE,fig.cap="波动率拟合图",cache=F,dev="png",fig.height=13,results='markup'----
knitr::include_graphics("./vol.png")

