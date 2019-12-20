## ----setup, echo=F-------------------------------------------------------
knitr::opts_knit$set(root.dir = getwd())
knitr::opts_chunk$set(echo = FALSE, results = 'hide')
knitr::opts_chunk$set(warning = FALSE, message=FALSE)


## ----prepare-------------------------------------------------------------
rm(list=ls())
options(digits=4)
options(scipen=100)
graphics.off()
Sys.setlocale("LC_ALL", "Chinese")


## ----vine,echo=FALSE,fig.cap="五维R藤模型结构图",cache=F,dev="png",results='markup'----
knitr::include_graphics("./vine.png")


## ----import-data---------------------------------------------------------
rm(list=ls())
## getwd()
## path='F:/github_paper/vine-copula-Zeng'
## setwd(path)
## load package
library(readxl)
library(xlsx)
library(splines)
library(fBasics)
library(tseries)
library(forecast)
library(fGarch)
library(nleqslv)
library(knitr)
library(FinTS)
library(rugarch)
library(VineCopula)
## load data
norisk_rate=read.csv('./no_risk_rate.csv',header=T)
no_risk_rate=norisk_rate[,3];date=as.Date(norisk_rate[,2])
stock_test=read.xlsx2('./sample_data.xlsx',sheetIndex=1,header=T,startRow = 1,endRow = 2359,
                      colClasses=c('character','character','numeric','numeric','numeric','character','character',rep('numeric',5)))
debt_book_value_test=read.xlsx2('./sample_data.xlsx',sheetIndex=2,header=T,startRow = 1,endRow = 42,
                      colClasses=c(rep('character',3),'numeric'))


## ----norisk-rate,eval=T, fig.cap="近十年一年期定期整存整取利率", dev='png'----
plot(date,no_risk_rate,type='l',ylim=c(0,4))


## ----cubicspline---------------------------------------------------------
n=length(stock_test$Stkcd)
  ## cubic spline
debt_book_value=spline(debt_book_value_test$Accper,debt_book_value_test$debt_book_value,n+1)
debt_book_value_date=as.Date(debt_book_value_test$Accper)


## ----cubspline,eval=T,fig.cap ="平安银行插值后日债务账面价值", dev='png'----
plot(debt_book_value_test$debt_book_value,type="p",xlab="日期",ylab="日债务账面价值",xaxt='n')
axis(1, at=1:41,labels = debt_book_value_test$Accper)
lines(debt_book_value$x,debt_book_value$y)


## ----return--------------------------------------------------------------
price=stock_test[,'price_new']
date=as.Date(stock_test[,'Trddt'])
price_rt=diff(log(price))


## ----log-return,eval=T,fig.cap = "平安银行近十年日收盘价和对数收益率时序图", dev='png'----
par(mfrow=c(2,1))
par(mar=c(2,4,1,2))
plot(date,price,type="l",main="(a)平安银行近十年收盘价时序图")
plot(date[-1],price_rt,type="l",cex.main=0.95,las=1,main='(b)平安银行近十年对数收益率时序图')


## ----price-rt-stats, results='markup'------------------------------------
price_rt_stats=basicStats(price_rt)
JB=jarque.bera.test(price_rt)
price_rt_stats_desc=as.data.frame(t(c(price_rt_stats[c(7,14,16,15),1],as.numeric(JB[[1]][1]))))
colnames(price_rt_stats_desc) <- c("均值","标准差","偏态系数","峰态系数","J-B统计量")

kable(price_rt_stats_desc,row.names =F,align = "c", caption="平安银行对数收益率描述性统计量",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----adf-test, results='markup'------------------------------------------
ADF=adf.test(price_rt)
PP=pp.test(price_rt)
adf_test=data.frame(a=c(ADF$statistic[1],PP$statistic[1]),b=c(0.0001,0.0001))
colnames(adf_test)=c("检验统计量的值","P值")
rownames(adf_test)=c("ADF检验","PP检验")

kable(adf_test,row.names =T,align = "c", caption="平安银行对数收益率平稳性检验",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----acf-pacf, fig.cap="Shibor收益率自相关图和偏自相关图", dev='png'-----
par(mfrow=c(2,1))
acf(price_rt,main="",xlab="滞后期",ylab="ACF",lag.max=20,ylim=c(-0.2,0.8))#画自相关图
title(main="(a)the ACF of Return",cex.main=0.95)
pacf(price_rt,main="",xlab="滞后期",ylab="PACF",lag.max=20,ylim=c(-0.2,0.2))#画偏自相关图
title(main="(b)the PACF of Return",cex.main=0.95)


## ----resi2, fig.cap="平安银行对数收益率去均值后残差平方时序图", dev='png'----
resi=price_rt-mean(price_rt)
resi2=(price_rt-mean(price_rt))^2
plot(date[-1],resi2,type='l',xlab='时间',ylab='对数收益率去均值后残差的平方')


## ----garch---------------------------------------------------------------
gjrgarch11_spec = ugarchspec(
  variance.model = list(model="gjrGARCH", garchOrder=c(1,1)),
  mean.model = list(armaOrder=c(0,0),include.mean=F))
garchfit11 = ugarchfit(spec=gjrgarch11_spec, data=price_rt)

gjrgarch12_spec = ugarchspec(
  variance.model = list(model="gjrGARCH", garchOrder=c(1,2)),
  mean.model = list(armaOrder=c(0,0),include.mean=F))
garchfit12 = ugarchfit(spec=gjrgarch12_spec, data=resi)

gjrgarch21_spec = ugarchspec(
  variance.model = list(model="gjrGARCH", garchOrder=c(2,1)),
  mean.model = list(armaOrder=c(0,0),include.mean=F))
garchfit21 = ugarchfit(spec=gjrgarch21_spec, data=resi)


## ----gjr-garch, results='markup'-----------------------------------------
garchfit11_coef=data.frame(garchfit11@fit$matcoef)
garchfit12_coef=data.frame(garchfit12@fit$matcoef)
garchfit21_coef=data.frame(garchfit21@fit$matcoef)
garchfit_coef=list(garchfit11_coef=garchfit11_coef,
                   garchfit12_coef=garchfit12_coef,
                   garchfit21_coef=garchfit21_coef)
for(j in 1:3){
  for(i in 1:nrow(garchfit_coef[[j]])){
    garchfit_coef[[j]][i,5]=paste(substr(garchfit_coef[[j]][i,1],1,7),
                           if(garchfit_coef[[j]][i,4]<=0.01){'***'}
                           else{if(garchfit_coef[[j]][i,4]>0.01 & garchfit_coef[[j]][i,4]<=0.05){'**'}
                             else{if(garchfit_coef[[j]][i,4]>0.05 & garchfit_coef[[j]][i,4]<=0.1){'*'}
                               else{''}}},
                           sep='',collapse = '')
  }
}
coef_result=data.frame(row.names=c('omega','alpha1','alpha2','beta1','beta2','gamma1','gamma2','AIC','BIC'))
coef_result[,1]=c(garchfit_coef$garchfit11_coef[1:2,5],'',
                  garchfit_coef$garchfit11_coef[3,5],'',
                  garchfit_coef$garchfit11_coef[4,5],'',signif(infocriteria(garchfit11)[1:2],5))
coef_result[,2]=c(garchfit_coef$garchfit12_coef[1:2,5],'',
                  garchfit_coef$garchfit12_coef[3:5,5],'',signif(infocriteria(garchfit12)[1:2],5))
coef_result[,3]=c(garchfit_coef$garchfit21_coef[1:4,5],'',
                  garchfit_coef$garchfit21_coef[5:6,5],signif(infocriteria(garchfit21)[1:2],5))
colnames(coef_result)=c('GJR-GARCH(1,1)','GJR-GARCH(1,2)','GJR-GARCH(2,1)')
kable(coef_result,row.names =T,align = "c", caption="平安银行对数收益率GJR-GARCH模型拟合结果",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----gjrgarch-resi-------------------------------------------------------
gjr_resi=residuals(garchfit11,standardize=T)
gjr_resi2=gjr_resi^2


## ----cca-solve-----------------------------------------------------------
### cca
##输入时间
r=stock_test$norisk_rate[-1]/100
T=1
B=debt_book_value$y[-c(1:2)]
##输入股权波动率和股权价值
EquityTheta=garchfit11@fit$sigma*sqrt(250)
E=stock_test$stock_market_value_new[-1]
##KMV 方程变形及求解
EtoB=E/B
x0=matrix(rep(1,2*nrow(stock_test)),nrow=nrow(stock_test))
z=list();Va=vector();AssetTheta=vector()
for (i in 1:(nrow(stock_test)-1)) {
  KMVfun=function(x){
    y=numeric(2);
    d1=(log(x[1]*EtoB[i])+(r[i]+0.5*x[2]^2)*T)/(x[2]*sqrt(T));
    d2=d1-x[2]*sqrt(T);
    y[1]=x[1]*pnorm(d1)-exp(-r[i]*T)*pnorm(d2)/EtoB[i]-1;
    y[2]=pnorm(d1)*x[1]*x[2]-EquityTheta[i];
    y
  }
  z[[i]]<-nleqslv(x0[i,], KMVfun, method="Newton") 
  Va[i]<-z[[i]]$x[1]*E[i]
  AssetTheta[i]<-z[[i]]$x[2]
}
##计算违约距离
DD=(Va-B)/(Va*AssetTheta)


## ----pinganbank-dd, fig.cap="平安银行近十年违约剧烈时序图", dev='png'----
datedd=as.Date(stock_test$Trddt[-1])
plot(datedd,DD,type='l')


## ----import-sector-dd----------------------------------------------------
data_dd=read.xlsx2("./sector_dd.xlsx",sheetIndex = 1,as.data.frame = TRUE,header = TRUE,
           colClasses = c("Date","numeric","numeric","numeric","numeric","numeric"))
dd_date=as.Date(data_dd[,1])
insurance=data_dd[,2]
multi_finance=data_dd[,3]
house=data_dd[,4]
stock=data_dd[,5]
bank=data_dd[,6]


## ----bank-sector-dd,eval=T,fig.cap = "银行板块近十年系统性风险时序图", dev='png'----
plot(dd_date,bank,type="l",ylim=c(0,6.5))


## ----multi-fin-sector-dd,eval=T,fig.cap = "多元金融板块近十年系统性风险时序图", dev='png'----
plot(dd_date,multi_finance,type="l",ylim=c(0,6.5))


## ----stock-sector-dd,eval=T,fig.cap = "券商信托板块近十年系统性风险时序图", dev='png'----
plot(dd_date,stock,type="l",ylim=c(0,6.5))


## ----house-sector-dd,eval=T,fig.cap = "房地产板块近十年系统性风险时序图", dev='png'----
plot(dd_date,house,type="l",ylim=c(0,6.5))


## ----insurance-sector-dd,eval=T,fig.cap = "保险板块近十年系统性风险时序图", dev='png'----
plot(dd_date,insurance,type="l",ylim=c(0,6.5))


## ----test-uniform-distribution-------------------------------------------
f_insurance=ecdf(data_dd[,2])
co_insurance=f_insurance(data_dd[,2])
f_multi_finance=ecdf(data_dd[,3])
co_multi_finance=f_multi_finance(data_dd[,3])
f_house=ecdf(data_dd[,4])
co_house=f_house(data_dd[,4])
f_stock=ecdf(data_dd[,5])
co_stock=f_stock(data_dd[,5])
f_bank=ecdf(data_dd[,6])
co_bank=f_bank(data_dd[,6])


## ----test-uniform, results='markup'--------------------------------------
test_uniform=data.frame(a=c(ks.test(co_insurance,'punif')[[1]],ks.test(co_insurance,'punif')[[2]]),
                    b=c(ks.test(co_multi_finance,'punif')[[1]],ks.test(co_multi_finance,'punif')[[2]]),
                    c=c(ks.test(co_house,'punif')[[1]],ks.test(co_house,'punif')[[2]]),
                    d=c(ks.test(co_stock,'punif')[[1]],ks.test(co_stock,'punif')[[2]]),
                    e=c(ks.test(co_bank,'punif')[[1]],ks.test(co_bank,'punif')[[2]]))
colnames(test_uniform)=c("保险板块","多元金融板块","房地产板块","券商信托板块","银行板块")
rownames(test_uniform)=c("K-S","P值")
kable(test_uniform,row.names =T,align = "c", caption="金融各板块系统性风险指标均匀分布检验",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----cor, results='markup'-----------------------------------------------
copula_data=data.frame(insurance=co_insurance,
                       multi_finance=co_multi_finance,
                       house=co_house,
                       stock=co_stock,
                       bank=co_bank)
cor0=cor(copula_data,method = 'kendall')
colnames(cor0)=c("保险板块","多元金融板块","房地产板块","券商信托板块","银行板块")
rownames(cor0)=c("保险板块","多元金融板块","房地产板块","券商信托板块","银行板块")
kable(cor0,row.names =T,align = "c", caption="金融各板块系统性风险指标相关系数矩阵",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----vine-first-tree,echo=FALSE,fig.cap="金融板块间系统性风险R藤模型树T1结构图",cache=F,dev="png",results='markup'----
knitr::include_graphics("./vine_first_tree.png")


## ----vine-first-tree-adj,echo=FALSE,fig.cap="调整后的金融板块间系统性风险R藤模型树T1结构图",cache=F,dev="png",results='markup'----
knitr::include_graphics("./vine_first_tree_adj.png")


## ----vine-fit------------------------------------------------------------
RVM=RVineStructureSelect(copula_data,familyset=c(0:5),progress=TRUE,rotations=FALSE,method='mle')
RVM2=RVineMLE(copula_data,RVM,maxit=200)


## ----vine-fit-fig,echo=FALSE,fig.cap="金融板块间系统性风险R藤模型结构图",cache=F,dev="png",results='markup'----
knitr::include_graphics("./vine_fit.png")


## ----rvm-summary1--------------------------------------------------------
RVM_summary1=data.frame(summary(RVM2$RVM)[,c(1:2,4:6)])


## ----vine-copula, results='markup'---------------------------------------
kable(RVM_summary1,row.names =F,align = "c", caption="板块间系统性风险最优Copula形式表",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----rvm-summary2--------------------------------------------------------
RVM_summary2=data.frame(summary(RVM2$RVM)[,c(1:2,8:9)])


## ----vine-ultd, results='markup'-----------------------------------------
kable(RVM_summary2,row.names =F,align = "c", caption="板块间系统性风险尾部相依系数表",
      longtable = TRUE, booktabs = TRUE, linesep="")


## ----rvm-summary3--------------------------------------------------------
RVM_kendall=data.frame(summary(RVM2$RVM)[,c(1:2,7)])


## ----vine-kendall, results='markup'--------------------------------------
kendall=c(cor0[4,3],cor0[4,2],cor0[4,1],cor0[1,5],cor0[3,2],cor0[3,1],cor0[4,5],cor0[2,1],cor0[3,5],cor0[2,5])
vine_kendall=cbind(RVM_kendall,kendall)
colnames(vine_kendall)=c("tree","edge","cond-tau","tau")
kable(vine_kendall,row.names =F,align = "c", caption="板块间系统性风险无条件Kendall相关系数与有条件Kendall相关系数对照表",
      longtable = TRUE, booktabs = TRUE, linesep="")

