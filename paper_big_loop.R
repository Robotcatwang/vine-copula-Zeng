rm(list=ls())
## load package
library(readxl)
library(splines)
library(fBasics)
library(tseries)
library(forecast)
library(fGarch)
library(nleqslv)
library(utils)
library(xlsx)
library(rJava)
## load data
getwd()
path='C:/Users/Administrator/Desktop'
setwd(path)

stock_test=read_excel('./test.xlsx',sheet=1,col_names=T,
                      col_types=c('text',"text",'numeric','numeric','numeric','text','text',rep('numeric',5)))
debt_book_value_test=read_excel('./test.xlsx',sheet=2,col_names=T,
                                col_types=c('text','text','text','numeric'))
delayday=read_excel('./test.xlsx',sheet=3,col_names=T,
                    col_types=c('text','numeric'))

stock_num=unique(stock_test$Stkcd)
stock_dd=list()
Pb=txtProgressBar(min=1,max=length(stock_num),style=3)
assign('last.warning',NULL,envir=baseenv())
warnings()

for(i in 244:length(stock_num)){
  options(warn=0)
  sample_stock=data.frame(stock_test[which(stock_test$Stkcd==stock_num[i]),])
  sample_stock_debt_book_value=data.frame(debt_book_value_test[which(debt_book_value_test$Stkcd==stock_num[i]),])
  n=nrow(sample_stock)
  zero=rep(0,n+delayday$day[i]+1);m0=vector();m1=vector()
  ## cubic spline
  debt_book_value=spline(c(1:length(sample_stock_debt_book_value[,'Accper'])),
                         sample_stock_debt_book_value[,'debt_book_value'],n+delayday$day[i]+1)
  if(length(which(debt_book_value$y<zero))>0){
    print(paste(i,'-',which(debt_book_value$y<zero)))
    for(k in 1:(n+delayday$day[i])){
      if(debt_book_value$y[k]>0 & debt_book_value$y[k+1]<=0){m0=c(m0,k)}
      if(debt_book_value$y[k]<=0 & debt_book_value$y[k+1]>0){m1=c(m1,k+1)}}
    for(g in 1:length(m0)){
      debt_book_value$y[(m0[g]+1):(m1[g]-1)]=
        mean(debt_book_value$y[m0[g]],debt_book_value$y[m1[g]])}
    print(paste(i,'-',which(debt_book_value$y<zero)))
    print(data.frame(m0=m0,m1=m1))
  }
  
  ## stock_volatility
  price=sample_stock[,'price_new']
  price_rt=diff(log(price))
  
  ### acf,pacf test
  fit=auto.arima(price_rt,max.p=3,max.q=3)
  
  ### garch
  gjrgarch11_spec = ugarchspec(
    variance.model = list(model="gjrGARCH", garchOrder=c(1,1)),
    mean.model = list(armaOrder=c(fit$arma[1],fit$arma[2]),include.mean=F))
  garchfit11 = ugarchfit(spec=gjrgarch11_spec, data=price_rt)
  
  gjrgarch12_spec = ugarchspec(
    variance.model = list(model="gjrGARCH", garchOrder=c(1,2)),
    mean.model = list(armaOrder=c(fit$arma[1],fit$arma[2]),include.mean=F))
  garchfit12 = ugarchfit(spec=gjrgarch12_spec, data=price_rt)
  
  gjrgarch21_spec = ugarchspec(
    variance.model = list(model="gjrGARCH", garchOrder=c(2,1)),
    mean.model = list(armaOrder=c(fit$arma[1],fit$arma[2]),include.mean=F))
  garchfit21 = ugarchfit(spec=gjrgarch21_spec, data=price_rt)
  
  if(mean(infocriteria(garchfit11)[1:2])<=mean(infocriteria(garchfit12)[1:2]) & 
     mean(infocriteria(garchfit11)[1:2])<=mean(infocriteria(garchfit21)[1:2])){
    stock_volatility=garchfit11@fit$sigma
  }else{
    if(mean(infocriteria(garchfit12)[1:2])<=mean(infocriteria(garchfit11)[1:2]) & 
       mean(infocriteria(garchfit12)[1:2])<=mean(infocriteria(garchfit21)[1:2])){
      stock_volatility=garchfit12@fit$sigma
    }else{
      stock_volatility=garchfit21@fit$sigma
    }
  }
  ## cca
  ## 输入时间
  r=sample_stock[,'norisk_rate'][-1]/100
  T=1
  B=debt_book_value[[2]][-c(1:(delayday$day[i]+2))]
  ##输入股权波动率和股权价值
  EquityTheta=stock_volatility*sqrt(250)
  E=sample_stock[,'stock_market_value_new'][-1]
  ##KMV 方程变形及求解
  EtoB=E/B
  m=n-1
  x0=matrix(rep(1,2*m),nrow=m)
  z=list();Va=vector();AssetTheta=vector()
  for (j in 1:m) {
    KMVfun=function(x){
      y=numeric(2);
      d1=(log(x[1]*EtoB[j])+(r[j]+0.5*x[2]^2)*T)/(x[2]*sqrt(T));
      d2=d1-x[2]*sqrt(T);
      y[1]=x[1]*pnorm(d1)-exp(-r[j]*T)*pnorm(d2)/EtoB[j]-1;
      y[2]=pnorm(d1)*x[1]*x[2]-EquityTheta[j];
      y
    }
    z[[j]]=nleqslv(x0[j,], KMVfun, method="Newton") 
    Va[j]=z[[j]]$x[1]*E[j]
    AssetTheta[j]=z[[j]]$x[2]
  }
  ##计算违约距离
  DD=(Va-B)/(Va*AssetTheta)
  stock_dd[[i]]=data.frame(stock_num=sample_stock[,'Stkcd'][-1],
                           stock_dt=sample_stock[,'Trddt'][-1],
                           assetvalue=Va,
                           DD=DD)
  #设置进度条
  Sys.sleep(0.02)
  setTxtProgressBar(Pb,i)
  if(length(warnings())>0){break}
}

jgc <- function(){
  gc()
  .jcall("java/lang/System", method = "gc")}
options(java.parameters = "-Xmx8000m")
stock_dd_new=stock_dd[[1]]
for(a in 2:length(stock_dd)){
  stock_dd_new=rbind(stock_dd_new,stock_dd[[a]])
}
write.table(stock_dd_new,'./stock_dd.csv',sep=',',col.names=TRUE,row.names=FALSE)
