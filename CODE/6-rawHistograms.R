### mTOR project


## load libraries
library(here)
library(caret)


### load data
getwd()
data<-read.csv(file=here("1-Data","MAK_rapamycin_CCM_20200530_nmol_g.csv"),header=TRUE)

###subset the data
colnames(data)
str(data, list.len=ncol(data))

list_factors<-c("Diet", "Number","Mouse.ID","Group")
data[,list_factors]<-lapply(data[,list_factors], as.factor)


colnames(data)
###create the histograms using a loop
for(i in 7:ncol(data))
{
  file.name<-paste0(colnames(data)[i],"_histogram.jpeg")
  jpeg(filename=file.name)
  hist(data[,i],xlab = colnames(data)[i], main = paste0("Histogram of ",colnames(data)[i]))
  dev.off()
}




###log transform data
idata.log<-log10(data[,-c(1:6)]+1)
idata<-cbind(data[,c(1:6)],idata.log)
colnames(idata)

###create the histograms using a loop
for(i in 7:ncol(data))
{
  file.name<-paste0(colnames(idata)[i],"_log_histogram.jpeg")
  jpeg(filename=file.name)
  hist(idata[,i],xlab = colnames(idata)[i], main = paste0("Histogram of Log ",colnames(idata)[i]))
  dev.off()
}




#Boxcox Transform
trans_setup<-caret::preProcess(data[,-c(1:6)], method=c("BoxCox"))
idata_trans<-predict(trans_setup,data[,-c(1:6)])

idata_trans<-cbind(data[,c(1:6)],idata_trans)
colnames(idata_trans)


###create the histograms using a loop
for(i in 7:ncol(data))
{
  file.name<-paste0(colnames(idata_trans)[i],"_boxCox_histogram.jpeg")
  jpeg(filename=file.name)
  hist(idata_trans[,i],xlab = colnames(idata_trans)[i], main = paste0("Histogram of BoxCox ",colnames(idata_trans)[i]))
  dev.off()
}



write.csv(data, file=here("1-Data","MAK_rapamycin_CCM_20200530_nmol_g_BoxCox.csv"))
