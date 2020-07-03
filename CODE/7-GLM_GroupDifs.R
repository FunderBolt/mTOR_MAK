### mTOR project


## load libraries
library(here)
#library(multcomp)

### load data
getwd()
data<-read.csv(file=here("1-Data","MAK_rapamycin_CCM_20200530_nmol_g_BoxCox.csv"),header=TRUE)

###subset the data
colnames(data)
str(data, list.len=ncol(data))

list_factors<-c("Diet", "Number","Mouse.ID","Group")
data[,list_factors]<-lapply(data[,list_factors], as.factor)


# # make a summary
# summary(data$Group)
# 
# 
# colnames(data)[1:20]
# fit<-glm(AMP_nmol_g~ Group, data=data, na.action = na.exclude)
# summary(fit)
# 
# 
# 
# 
# fit_mcomp<-glht(fit, mcp(Group="Tukey"))
# names(fit_mcomp)
# 
# fit_mcomp$coef
# out_comp<-summary(glht(fit, mcp(Group="Tukey")))
# names(out_comp)
# 
# cbind(out_comp$test$coefficients, out_comp$test$sigma, out_comp$test$pvalues)

################################################## 3 groups - 5 samples
#####################################
colnames(data)[1:20]

idata<-data[!data$Group==3,]
idata<-data[!data$Group==2,]
idata<-data[!data$Group==1,]
###GLMs for all metabolites --- create a function
GLM.run<-function(y) {
  form <- as.formula(paste0(y,"~ Group"))
  fit<-(glm(form, data=idata, na.action = na.exclude))
}

### apply the function to all metabolites
GLMs.out <- lapply(colnames(idata)[c(8:ncol(idata))],GLM.run )

#### print out results
results<-lapply(GLMs.out, function(x){summary(x)})
results


#Pull coefficients from models
estim.coef.results<-sapply(results, function(x){coef(x)})
estim.coef.results

## transpose the table
estim.coef.results<-t(estim.coef.results)
estim.coef.results


#the output of above commands create a matrix, but dataframes are better to indicate column names and rownames
class(estim.coef.results)
estim.coef.results<-as.data.frame(estim.coef.results)

#add rownames for identifiers
rownames(estim.coef.results)<-colnames(idata)[c(8:ncol(idata))]
head(estim.coef.results)

#add columns names
results[1]
colnames(estim.coef.results)<-c("est.intercept","est.G1_G2", 
                                "se.intercept","se.G1_G2",
                                "t.intercept","t.G1_G2",
                                "p.intercept","p.G1_G2")
head(estim.coef.results)


colnames(estim.coef.results)<-c("est.intercept","est.G1_G3", 
                                "se.intercept","se.G1_G3",
                                "t.intercept","t.G1_G3",
                                "p.intercept","p.G1_G3")
head(estim.coef.results)


colnames(estim.coef.results)<-c("est.intercept","est.G2_G3", 
                                "se.intercept","se.G2_G3",
                                "t.intercept","t.G2_G3",
                                "p.intercept","p.G2_G3")
head(estim.coef.results)
#ALL_estim.coef.results<-estim.coef.results
ALL_estim.coef.results<-cbind(ALL_estim.coef.results,  estim.coef.results)


ALL_estim.coef.results



### add in FDR corrected p-values
head(ALL_estim.coef.results)
p.vals<-c(ALL_estim.coef.results$p.G1_G2, ALL_estim.coef.results$p.G1_G3, ALL_estim.coef.results$p.G2_G3)


# apply correction -- > some correction option methods: "bonferroni", "fdr"
final<-p.adjust(p.vals, method="fdr")
head(final)
final<-signif(final, digits=2)
final

#### add in FDR p-values
head(estim.coef.results)
ALL_estim.coef.results$FDR.G1_G2<-final[1:51]
ALL_estim.coef.results$FDR.G1_G3<-final[52:102]
ALL_estim.coef.results$FDR.G2_G3<-final[103:153]

head(ALL_estim.coef.results)
###


ALL_estim.coef.results<-apply(ALL_estim.coef.results,2, signif, digits=2)
head(ALL_estim.coef.results)

write.csv(ALL_estim.coef.results,file="BoxCoxTransform_coefficients_allCCs_round2digits.csv")








# ###Residual Plots
# Residual_plot.run <-function(y) {
#   form <- as.formula(paste0(y,"~ case_SB_CC + chrosam_wealth1 + chrosam_age_m + sex + chrosam_hiv_status_0negative_1positive_2unknown "))
#   fit<-glm(form, data=idata, na.action = na.exclude)
#   file.name<-paste0(y,".GLM_logTransform_residual_plot.jpeg")
#   jpeg(filename=file.name)
#   print(plot(fit,which=1))
#   dev.off()
# }
# 
# 
# setwd("C:\\Users\\Celine\\Dropbox\\Bandsma.Lab\\1.Projects\\2016_ChroSAM\\")
# lapply(colnames(idata)[15:228],Residual_plot.run )
# 
# ###QQ-plots
# qqplots.run <-function(y) {
#   form <- as.formula(paste0(y,"~ case_SB_CC + chrosam_wealth1 + chrosam_age_m + sex + chrosam_hiv_status_0negative_1positive_2unknown"))
#   fit<-glm(form, data=idata, na.action = na.exclude)
#   file.name<-paste0(y,".GLM_logTransform_qqplot.jpeg")
#   jpeg(filename=file.name)
#   print(plot(fit,which=2))
#   dev.off()
# }
# setwd("C:\\Users\\Celine\\Desktop\\ChroSAM\\GLM_qqplots")
# lapply(colnames(idata)[c(11:226)],qqplots.run )
# 
# 
# 


# ###Boxplots
# colnames(idata)[1:20]
# 
# #set location to where you want the plots saved
# setwd("C:\\Users\\Celine\\Desktop\\ChroSAM\\Boxplots")
# 
# #create all boxplots saved as jpeg
# for (i in 11:226)
# {
#   plot<-boxplot(idata[,i]~idata[,2], main= paste0(colnames(idata)[i]," Boxplot"), xlab = colnames(idata)[i])
#   file.name<-paste0(colnames(idata)[i],"_BoxPlots",".jpeg")
#   dev.copy(jpeg, file.name)
#   dev.off()
# }
# 
# 
# #create all boxplots saved as .svg files
# for (i in 11:226)
# {
#   plot<-boxplot(idata[,i]~idata[,5], main= paste0(colnames(idata)[i]," Boxplot"), xlab = colnames(idata)[i])
#   file.name<-paste0(colnames(idata)[i],"_BoxPlots",".svg")
#   dev.copy(svg, file.name)
#   dev.off()
# }
# 

