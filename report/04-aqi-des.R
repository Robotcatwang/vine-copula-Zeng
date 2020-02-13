## ----setup, echo=F, message=F--------------------------------------------

################# 第 4 章 R 程序代码  ####################


knitr::opts_knit$set(root.dir = getwd())
knitr::opts_chunk$set(echo = FALSE, results = 'hide')
knitr::opts_chunk$set(warning = FALSE, message=FALSE)
options(knitr.kable.NA = '')


## ----prepare-------------------------------------------------------------
rm(list=ls())
options(digits=4)
options(scipen=100)
graphics.off()
Sys.setlocale("LC_ALL", "Chinese")
windowsFonts(msyh=windowsFont("微软雅黑"))
library("kableExtra")


## ----VaR8, fig.cap="1%显著性水平下 GAS-t 回测图", dev='png'--------------
n.ots    <- 750
n.its    <- 1733 
alpha    <- 0.01 
k.update <- 100

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

dat=read.csv("./data/jiesuan.csv",header=T)
shibor<-dat[,2]
shibor<-shibor[1:2484]
shibor.date<-as.Date(dat[,1])
shibor.date<-shibor.date[1:2484]

#计算对数收益率
shibor.rt<-diff(log(shibor))
shibor.rt<-shibor.rt[1:2483]
shibor.rt = na.omit(shibor.rt)
fit=stats::arima(shibor.rt,order = c(1,0,1))
shibor.rt.r<-stats::residuals(fit)

InSampleData = shibor.rt.r[1:1733]
OutSampleData = shibor.rt.r[1734:2483]


herspec=UniGASSpec(Dist="std",ScalingType="Identity",GASPar=list(location=TRUE,scale=TRUE,shape=TRUE))

herfit=UniGASFit(herspec,data=InSampleData,fn.optimizer = fn.optim, Compute.SE = TRUE)

models=list(herspec)

var.res <- read.csv("./VaR.csv")
VaR_N <- var.res[,2]
ES_N <- var.res[,3]
y.ots <- as.vector(var.res[,4])

library("zoo")
time.index <- zoo::index(shibor.rt.r)[(n.its + 1):(n.ots + n.its)]
y_ots <- zoo::zoo(y.ots, order.by = time.index)
VaR   <- zoo::zoo(VaR_N, order.by = time.index)
ES   <- zoo::zoo(ES_N, order.by = time.index)

par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, xlab = "",
     ylab = "", col = "black", cex.axis = 1, cex.lab = 1,cex=1, pch = 19)
lines(VaR[,1], type = 'l', col = "green", lwd = 2, lty = "dashed")
lines(ES[,1], type = 'l', col = "red", lwd = 2)
legend("topleft", legend =c("VaR 1% - VaR", "VaR 1% -ES"), 
       col = c("green","red"), lwd = 2, cex = 1, lty = c("dashed","solid"))
abline(h = 0)#在y=0添加次要刻度线


## ----VaR9, fig.cap="5%显著性水平下 GAS-t 回测图", dev='png'--------------
var.res1 <- read.csv("./VaR1.csv")
VaR_N2 <- var.res1[,2]
ES_N2 <- var.res1[,3]
y.ots <- as.vector(var.res1[,4])

VaR2   <- zoo::zoo(VaR_N2, order.by = time.index)
ES2   <- zoo::zoo(ES_N2, order.by = time.index)

par(mfrow = c(1, 1))
plot(y_ots, type = 'p', las = 1, xlab = "",
     ylab = "", col = "black", cex.axis = 1, cex.lab = 1,cex=1, pch = 19)
lines(VaR2[,1], type = 'l', col = "green", lwd = 2, lty = "dashed")
lines(ES2[,1], type = 'l', col = "red", lwd = 2)
legend("topleft", legend =c("VaR 5% - VaR", "VaR 5% -ES"), 
       col = c("green","red"), lwd = 2, cex = 1, lty = c("dashed","solid"))
abline(h = 0)#在y=0添加次要刻度线


## ----egarch-var,fig.cap="1%显著性水平下 EGARCH-t 回测图",echo=FALSE,cache=F,dev="png",results='markup'----
knitr::include_graphics("./garchvar.png")


## ----es-forcast, results='markup'----------------------------------------
es_forcast=data.frame(a=c(-0.018278773,-0.017602375,-0.016910876,-0.016272798,-0.017580637,-0.018146741,-0.017741588,-0.016913585,-0.016216954,-0.016191264,-0.01810264,-0.017761975,-0.017105657),
                    b=c(-0.022741092,-0.021889395,-0.021036049,-0.020240589,-0.021796041,-0.022472494,-0.02198054,-0.020988671,-0.020157886,-0.020164734,-0.022445645,-0.022076299,-0.021269631),
                    c=c(-0.0049963,-0.007968853,0.008171182,-0.001285402,0.00176795,0.009809547,0.007761012,0.001458591,-0.004139269,-0.010945445,0.009625589,-0.005607119,-0.014411272))
colnames(es_forcast)=c("VaR","ES","实际值")
rownames(es_forcast)=c("1","2","3","4","5","6","7","8","9","10","11","12","13")
kable(es_forcast,row.names =T,align = "c", caption="GAS-t的VaR，ES及实际值",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----test8-tab,eval=T,echo=F---------------------------------------------
library("GAS")
library(rugarch)
library(fGarch)
library(parallel)
spec1=ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                 mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                 distribution.model = "std")
models1=list(spec1)

POF<-UC<-CC<-DQ<-UC.pval<-CC.pval <- DQ.pval <- vector("double", length(models))
var.res <- read.csv("./VaR.csv")
VaR <- var.res[,2:3]
y.ots <- as.vector(var.res[,4])

pof=function(x){
  k=0
  for(i in 1:750){
    if(x[i]>shibor.rt.r[1733+i]){k=k+1}
  }
  k/750
}

for (j in 1:length(models)) {
    test0 <- GAS::BacktestVaR(data  = y.ots,
                             VaR   = VaR[,j],
                             alpha = 0.01)
    UC[j] <- test0$LRuc[1]  
    CC[j] <- test0$LRcc[1] 
    DQ[j] <- test0$DQ$stat
    
    UC.pval[j] <- test0$LRuc[2]  
    CC.pval[j] <- test0$LRcc[2] 
    DQ.pval[j] <- test0$DQ$pvalue
    
    POF[j]=pof(VaR[,j])} 
test0.print=t(data.frame(POF,UC,UC.pval,CC,CC.pval,DQ,DQ.pval))
test0.print <- cbind(c("失败率","Kupiec检验统计量","Kupiec检验P值","条件覆盖测试统计量",
                        "条件覆盖测试P值","动态分位数测试统计量","动态分位数测试P值"), test0.print)
colnames(test0.print)=c("统计量","VaR")
rownames(test0.print)=c("失败率","Kupiec检验统计量","Kupiec检验P值","条件覆盖测试统计量",
                        "条件覆盖测试P值","动态分位数测试统计量","动态分位数测试P值")

write.csv(test0.print, "test8.csv", row.names = F)




## ----test8, eval=T,results='markup', cache=F-----------------------------
tab8 <- read.csv('./test8.csv')
knitr::kable(tab8, row.names =F, align = "l", caption="GAS-t回测检验",
      longtable = TRUE, booktabs = TRUE, linesep  = "", escape = F)


## ----test9-tab,eval=T,echo=F---------------------------------------------
library("GAS")
library(rugarch)
library(fGarch)
library(parallel)
spec1=ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                 mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                 distribution.model = "std")
models1=list(spec1)
POF<-UC<-CC<-DQ<-UC.pval<-CC.pval <- DQ.pval <- vector("double", length(models1))
var.res <- read.csv("./VaR2.csv")
VaR <- var.res[,2:3]
y.ots <- as.vector(var.res[,3])

pof=function(x){
  k=0
  for(i in 1:750){
    if(x[i]>shibor.rt.r[1733+i]){k=k+1}
  }
  k/750
}

for (j in 1:length(models1)) {
    test0 <- GAS::BacktestVaR(data  = y.ots,
                             VaR   = VaR[,j],
                             alpha = 0.01)
    
    UC.pval[j] <- test0$LRuc[2]  
    CC.pval[j] <- test0$LRcc[2] 
    DQ.pval[j] <- test0$DQ$pvalue
    
    POF[j]=pof(VaR[,j])} 
test0.print=t(data.frame(POF,UC.pval,CC.pval,DQ.pval))
test0.print <- cbind(c("失败率","Kupiec检验P值","条件覆盖测试P值", "动态分位数测试P值"), test0.print)
colnames(test0.print)=c("统计量", "VaR")
rownames(test0.print)=c("失败率","Kupiec检验P值","条件覆盖测试P值", "动态分位数测试P值")

write.csv(test0.print, "test9.csv", row.names = F)




## ----test9, eval=T,results='markup', cache=F-----------------------------
tab9 <- read.csv('./test9.csv')
knitr::kable(tab9, row.names =F, align = "l", caption="EGARCH-t回测检验",
      longtable = TRUE, booktabs = TRUE, linesep  = "", escape = F)


## ----forcast-var, results='markup'---------------------------------------
forcast_var=data.frame(a=c(1.34E-05,2.04E-05,3.01E-06,4.34E-06,0.000201204)*1E+03,
                    b=c(1.34E-05,1.52E-05,1.91E-05,2.36E-05,2.92E-05)*1E+03,
                    c=c(0.00E+00,5.20E-06,-1.61E-05,-1.93E-05,1.72E-04)*1E+03,
                    d=c(-0.024882357,-0.026127331,-0.02655962,-0.027098548,-0.028248224),
                    e=c(-0.030691363,-0.032252327,-0.032800654,-0.033595117,-0.035159475))
colnames(forcast_var)=c("波动率真实值","波动率预测值","波动误差","VaR","ES")
rownames(forcast_var)=c("2485(2019.03.25)","2486(2019.03.26)","2487(2019.03.27)","2488(2019.03.28)","2489(2019.03.29)")
kable(forcast_var,row.names =T,align = "c", caption="未来五天黄金期货价格预测情况",
      longtable = TRUE, booktabs = TRUE, linesep="")

