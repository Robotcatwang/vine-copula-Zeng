## ----setup, echo=F, message=F--------------------------------------------

################# 第 4 章 R 程序代码  ####################


knitr::opts_knit$set(root.dir = getwd())
knitr::opts_chunk$set(echo = FALSE, results = 'hide')
knitr::opts_chunk$set(warning = FALSE, message=FALSE)
options(knitr.kable.NA = '')

### --prepare----------------------------------------
rm(list=ls())
options(digits=4)
options(scipen=100)
graphics.off()
Sys.setlocale("LC_ALL", "Chinese")
windowsFonts(msyh=windowsFont("微软雅黑"))
library("kableExtra")


### --import-data2-----------------------------
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

### --import-sector-dd2----------------------------------
data_dd=read.xlsx2("./sector_dd.xlsx",sheetIndex = 1,as.data.frame = TRUE,header = TRUE,
                   colClasses = c("Date","numeric","numeric","numeric","numeric","numeric"))
dd_date=as.Date(data_dd[,1])
insurance=data_dd[,2]
multi_finance=data_dd[,3]
house=data_dd[,4]
stock=data_dd[,5]
bank=data_dd[,6]


### --test-uniform-distribution--------------------------
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


### --test-uniform, results='markup'-------------------
test_uniform=data.frame(a=c(ks.test(co_insurance,'punif')[[1]],ks.test(co_insurance,'punif')[[2]]),
                        b=c(ks.test(co_multi_finance,'punif')[[1]],ks.test(co_multi_finance,'punif')[[2]]),
                        c=c(ks.test(co_house,'punif')[[1]],ks.test(co_house,'punif')[[2]]),
                        d=c(ks.test(co_stock,'punif')[[1]],ks.test(co_stock,'punif')[[2]]),
                        e=c(ks.test(co_bank,'punif')[[1]],ks.test(co_bank,'punif')[[2]]))
colnames(test_uniform)=c("保险板块","多元金融板块","房地产板块","券商信托板块","银行板块")
rownames(test_uniform)=c("K-S","P值")
kable(test_uniform,row.names =T,align = "c", caption="金融各板块系统性风险指标均匀分布检验",
      longtable = TRUE, booktabs = TRUE, linesep="")


### --cor, results='markup'--------------------------
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


### --vine-first-tree,echo=FALSE,fig.cap="金融板块间系统性风险R藤模型树T1结构图"------------------
knitr::include_graphics("./vine_first_tree.png")


### --vine-first-tree-adj,echo=FALSE,fig.cap="调整后的金融板块间系统性风险R藤模型树T1结构图"------------------
knitr::include_graphics("./vine_first_tree_adj.png")


### --vine-fit------------------------
RVM=RVineStructureSelect(copula_data,familyset=c(0:5),progress=TRUE,rotations=FALSE,method='mle')
RVM2=RVineMLE(copula_data,RVM,maxit=200)


### --vine-marix-----------------------------------
vine_matrix=data.frame(RVM2$RVM$Matrix) 
rownames(vine_matrix)=c("A","B","C","D","E") 
colnames(vine_matrix)=c("A","B","C","D","E") 
kable(vine_matrix,row.names =T,align = "c", caption="各板块系统性风险藤结构矩阵表达式",    longtable = TRUE, booktabs = TRUE, linesep="") 


### --vine-fit-fig,echo=FALSE,fig.cap="金融板块间系统性风险R藤模型结构图"--------------------------
knitr::include_graphics("./vine_fit.png")


### --rvm-summary1-----------------------------------
RVM_summary1=data.frame(summary(RVM2$RVM)[,c(1:2,4:6)])


### --vine-copula-----------------------------
kable(RVM_summary1,row.names =F,align = "c", caption="板块间系统性风险最优Copula形式表",
      longtable = TRUE, booktabs = TRUE, linesep="")


### --rvm-summary2--------------------------
RVM_summary2=data.frame(summary(RVM2$RVM)[,c(1:2,8:9)])


### --vine-ultd----------------------------------
kable(RVM_summary2,row.names =F,align = "c", caption="板块间系统性风险尾部相依系数表",
      longtable = TRUE, booktabs = TRUE, linesep="")


### --rvm-summary3------------------------------
RVM_kendall=data.frame(summary(RVM2$RVM)[,c(1:2,7)])


### --vine-kendall--------------------------
kendall=c(cor0[4,3],cor0[4,2],cor0[4,1],cor0[1,5],cor0[3,2],cor0[3,1],cor0[4,5],cor0[2,1],cor0[3,5],cor0[2,5])
vine_kendall=cbind(RVM_kendall,kendall)
colnames(vine_kendall)=c("tree","edge","cond-tau","tau")
kable(vine_kendall,row.names =F,align = "c", caption="板块间系统性风险无条件Kendall相关系数与有条件Kendall相关系数对照表",
      longtable = TRUE, booktabs = TRUE, linesep="")
