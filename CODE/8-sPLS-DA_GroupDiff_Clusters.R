### MixOmics: sPLS-DA script
### mTOR project


## load libraries
library(mixOmics)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(here)
library(mgsub)
library(ggpubr)
library(ggforce)
library(ggROC)
library(pROC)
library(directlabels)


### load data
getwd()
data<-read.csv(file=here("1-Data","MAK_rapamycin_CCM_20200530_nmol_g_BoxCox.csv"),header=TRUE, row.names=1)
data<-read.csv(file=here("1-Data","MAK_rapamycin_CCM_20200530_protcorr_BoxCox.csv"),header=TRUE, row.names=1)

###subset the data
colnames(data)
str(data, list.len=ncol(data))

list_factors<-c("Diet", "Number","Mouse.ID","Group")
data[,list_factors]<-lapply(data[,list_factors], as.factor)
################### check coding of variables
str(data, list.len=ncol(data))



###############################
## running PCA

#pull data of interest
colnames(data)

### all variables
PCA1<-PCA(data[,c(4,7:ncol(data))],quali.sup =c(1), graph=T)
PCA1

plot(PCA1, habillage=1, col.hab=c("red","black","green"), label="none")






#####################################################################
#####################################################################
### Running PLS-DA

head(data)
X<-as.matrix(data[,c(7:ncol(data))])

### Transform all variables --> no need loaded boxcox transformed
#X<-log(X+10)


#set outcome
outcome<-data$Group


### plsda between groups
out.plsDA<-plsda(X,outcome, ncomp=2, scale=TRUE, mode="regression")
out.plsDA$explained_variance

perf.out<-mixOmics::perf(out.plsDA, validation = c("loo"), progressBar = TRUE, auc=TRUE)
summary(perf.out)

 #get.confusion_matrix(truth = data$Group, predicted=perf.out$predict)
 perf.out$error.rate
 perf.out$error.rate.class
 perf.out$predict
 perf.out$class
 plot(perf.out, overlay = 'measure', sd=TRUE)
 #plot(perf.out$error.rate.class, overlay = 'measure', sd=TRUE)

 plotIndiv(out.plsDA,cex=5)
 plotVar(out.plsDA,cex=5)
 
 
 ################################################################
 # create list of variables to keep and test
 list.keepX<-c(seq(2,30,2))
 
 ###### tune splsda
 colnames(data)[1:20]
 out.splsDA<-mixOmics::tune.splsda(X,outcome, test.keepX=list.keepX, ncomp=2, validation = "loo", max.iter = 500, tol = 1e-06,
                                   dist = "centroids.dist", measure = "BER")
 
 out.splsDA$choice.keepX
 plot(out.splsDA)
 
 
 ######
 ### run plsda 
 colnames(idata)[1:20]
 out.splsDA<-mixOmics::splsda(X,outcome, keepX=c(20,20), ncomp=2, max.iter = 500, tol = 1e-06)
 out.splsDA$explained_variance
 
 ### get performance values
 perf.out<-perf(out.splsDA, validation=c("loo"), criterion="all", auc=T)
 perf.out
 
 plot(perf.out)
 perf.out$error.rate
 perf.out$error.rate.sd
 
 perf.out$error.rate.class
 
 auroc(out.plsDA, roc.comp = 1)
 auroc(out.plsDA, roc.comp = 2)
 
 ##AUC
 #perf.out$auc # average of repetitions
 #perf.out$auc.all # AUC for each repeat
 
 
 
 plotIndiv(out.splsDA,cex=5)
 
 ###########################################
 #### plotting variables
 plotVar(out.splsDA,cex=3)
 plotVar(out.splsDA,cex=3, cutoff=0.6)
 
 
 ### plotting stability of variables
 nrow(perf.out$features$stable[[1]])
 plot(perf.out$features$stable[[1]][1: nrow(perf.out$features$stable[[1]])], type = 'h', ylab = 'Stability', 
      xlab = 'Features', main = 'Comp 1', las =2, cex=0.3)
 
 
 nrow(perf.out$features$stable[[2]])
 plot(perf.out$features$stable[[2]][1: nrow(perf.out$features$stable[[1]])], type = 'h', ylab = 'Stability', 
      xlab = 'Features', main = 'Comp 2', las =2, cex=0.3)
 
 
 
 
 # collect loadings into file
 
 feature_loadings<-as.data.frame(out.plsDA$loadings$X)
 feature_loadings$Metabolite<-rownames(feature_loadings)
 colnames(feature_loadings)<-c("PLS_Loading_comp1","PLS_Loading_comp2","Metabolite")
 head(feature_loadings)
 
 
 
 feature_loadings_sPLS<-as.data.frame(out.splsDA$loadings$X)
 feature_loadings_sPLS$Metabolite<-rownames(feature_loadings_sPLS)
 colnames(feature_loadings_sPLS)<-c("sPLS_Loading_comp1","sPLS_Loading_comp2","Metabolite")
 head(feature_loadings_sPLS)

 
 plsda_features<-left_join(feature_loadings, feature_loadings_sPLS, by="Metabolite")
 
 
 variables.vip<-data.frame(vip(out.plsDA))
 variables.vip$Metabolite<-rownames(variables.vip)
 colnames(variables.vip)<-c("PLS_VIP_comp1","PLS_VIP_comp2","Metabolite")
 head(variables.vip)
 
 
 plsda_features<-left_join(plsda_features, variables.vip, by="Metabolite")
 head(plsda_features)
 
 
 
 res.perf<-as.data.frame(perf.out$features$stable$comp1)
 head(res.perf)
 colnames(res.perf)<-c("Metabolite","Freq_comp1")
 
 
 plsda_features<-left_join(plsda_features, res.perf, by="Metabolite")
 
 
 res.perf<-as.data.frame(perf.out$features$stable$comp2)
 head(res.perf)
 colnames(res.perf)<-c("Metabolite","Freq_comp2")
 
 

 plsda_features<-left_join(plsda_features, res.perf, by="Metabolite")
 head(plsda_features)
 
 


 variables.vip<-data.frame(vip(out.splsDA))
 variables.vip$Metabolite<-rownames(variables.vip)
 colnames(variables.vip)<-c("sPLS_VIP_comp1","sPLS_VIP_comp2","Metabolite")
 head(variables.vip)
 
 
 plsda_features<-left_join(plsda_features, variables.vip, by="Metabolite")
 head(plsda_features)
 
 

 
 write.csv(plsda_features,file=paste0("plsDA_features_Comp1&2_", Sys.Date() , ".csv"))
 write.csv(plsda_features,file=paste0("plsDA_features_Comp1&2_protcorr_", Sys.Date() , ".csv"))
 
 

 
############################################# Plotting results

### plotting vars with 0.6 cut off
vars.in<-plotVar(out.splsDA, plot=TRUE,cex=3,cutoff=0.5)
vars.in

#### list to include
metabs<-plsda_features[ plsda_features$sPLS_VIP_comp1 > 1 | plsda_features$sPLS_VIP_comp2 > 1,  ]



plot_Data<-vars.in[,c("names","x","y")]
### restict data to those that meet cut off
plot_Data<-plot_Data[plot_Data$names %in% metabs$Metabolite, ]

## cleaning variable names 
dput(plot_Data$names)
plot_Data$names<-c("GMP", "GDP", "cyclic.AMP", "Acetyl Phosphate", 
                   "2,3 bisP Glycerate", "Sedoheptulose 7P", "Acetyl CoA", "NAD", 
                   "NADH", "NADPH", "Erythrose 4P", "Glycerylaldehyde 3P", "Mannose 6P", 
                   "Ribose 5P", "a-Hydroxyglutaric ac.", "Glycolic ac.", "Isocitric ac.")


plot_Data$names<-c("GDP", "UDP", "cyclic.AMP", "Acetyl.Phosphate", "X2.3.bisP.Glycerate", 
  "total.of.Fructose.bisP.Glucose.1.6.bisP", "Sedoheptulose.7P", 
  "Acetyl.CoA", "Hs.CoA", "NAD.", "NADH", "NADPH", "Glucose", "Glycerylaldehyde.3P", 
  "Ribose.5P", "a.Hydroxyglutaric.acid", "Citric.acid", "Glycolic.acid", 
  "Isocitric.acid", "Pyruvic.acid", "Succinic.acid")

colnames(plot_Data)



p1<-ggplot(aes(x=x, y=y), data=plot_Data)+
  geom_point()+
 scale_x_continuous(lim=c(-1.5,1.5))+
  scale_y_continuous(lim=c(-1.5,1.5))+
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), inherit.aes = FALSE, colour="grey30", size=0.75)+
  geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), inherit.aes = FALSE, colour="grey50", size=0.5, linetype=2)+
  geom_hline(yintercept=0, linetype="dashed",color="grey") +
  geom_vline(xintercept=0, linetype="dashed", color="grey") +
  geom_segment(aes(xend = x, yend = y, x = 0, y = 0),size = 0.75, colour = "black", data = plot_Data, 
               arrow = arrow(length = unit(10, "points"), type = "open", angle = 25))+
  geom_dl(aes(x = x, y = y, label = names), 
          method = list(dl.trans(x = x - 0.2, y = y + 0), "first.bumpup", cex = 0.5)) +
  labs(x = "Component-1", y = "Component-2")+
  theme_bw()+
  theme(axis.title = element_text(size = 14), axis.text=element_text(color = "grey30", size=14))
p1




vars.in<-plotIndiv(out.plsDA)[[1]]


p2<-ggplot(aes(x=x, y=y), data=vars.in)+
  geom_point(aes(col=group, fill=group), shape=c(21), size=6)+
  scale_shape_manual(values=c(19,21,3,19,21,3))+
  #scale_x_continuous(lim=c(0,2200))+
  #scale_y_continuous(lim=c(0,130))+
  scale_color_manual(values=c("black","black","grey60"))+
  scale_fill_manual(values=alpha(c("white","black","grey60"), 1))+
  labs(x = paste0("Variate-1 (", 100*round(out.splsDA$explained_variance$X[1],2) , "% explained variance)"), y = paste0("Variate-2 (", 100*round(out.splsDA$explained_variance$X[2],2) , "% explained variance)"))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), axis.text=element_text(color = "grey30", size=14), legend.position = "right")+
  stat_ellipse(aes(group = na.omit(group)), na.rm = FALSE, show.legend = NA, level = 0.95) +
  #stat_conf_ellipse(aes(group = na.omit(group)), na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, level = 0.95, npoint = 100, bary = TRUE) +
  stat_mean(aes(group = group),size=5, shape=8, color= c("black","black","grey30"), na.rm = TRUE, show.legend = NA,
            inherit.aes = TRUE)
p2


p_all<-ggarrange(p1,p2, ncol = 2, nrow =1, common.legend = FALSE) 
p_all

# saves as .svg size 1000 800
ggsave(file="Fig_PLSDA_byGroups.svg", plot=p_all, width=12, height=6)
ggsave(file="Fig_PLSDA_byGroups.pdf", plot=p_all, width=12, height=6)



# saves as .svg size 1000 800
ggsave(file="Fig_PLSDA_byGroups_noProtcorr.svg", plot=p_all, width=12, height=6)
ggsave(file="Fig_PLSDA_byGroups_noProtcorr.pdf", plot=p_all, width=12, height=6)


# #### Making Roc curves
# colnames(data)
# 
# #####################################################################
# 
# auc_groups<-auroc(out.plsDA,roc.comp=1)
# auc_groups_data<-auc_groups$graph.Comp1$data
# 
# 
# #####################################################################
# auc_Groups_1_2_3<-ggplot()+
#   geom_line(aes(x=Specificity, y=Sensitivity), color="grey10",size = 1.25, alpha=0.75, data=auc_groups_data)+
#   geom_line(aes(x=Specificity, y=Sensitivity),color="grey60", size = 1.25, alpha=0.75, data=auc_groups_data)+
#   geom_segment(aes(xend = 100, yend = 100, x = 0, y = 0), size = 0.75, colour = "grey60", linetype="dashed")+
#   labs(title= "sW vs. CP", x = "100 - Specificity (%)", 
#        y = "Sensitivity (%)")+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         title = element_text(size = 16),
#         axis.title = element_text(size = 14), 
#         axis.text=element_text(color = "grey30", size=14), legend.position = "right")
# 
# auc_Groups_1_2_3
# 
# 
# 
# # saves as .svg size 1000 800
# ggsave(file="Fig3_PLSDA_anthroOnly_vs_withBIVA_AUC.svg", plot=p_all_auc, width=10, height=5)
# ggsave(file="Fig3_PLSDA_anthroOnly_vs_withBIVA_AUC.pdf", plot=p_all_auc, width=10, height=5)
# 
# 
