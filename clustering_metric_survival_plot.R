library(ggplot2)
library(ggpubr)

library(pROC)
library(survival)
library(survminer)
library(fpc)
par(mfrow=c(2,2))
setwd("./Analysis_brca/icluster/")

#Clustring Metric (Calinski Harabasz criteria  and silhouette index) Visualization

#Clustric metric visualization for sclae.lambda=c(1,1,1)
output=alist()
files=grep("cv.ics_1_1_1",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda) 
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)

setwd("C:/Users/pavel/Documents/")
clin_in_new<-read.table("./Analysis_brca/data/clin_brca_use.txt")
clin_in_new$deceased=clin_in_new$vital_status=="Dead"
clin_in_new$overall_survival=ifelse(clin_in_new$deceased,clin_in_new$days_to_death, clin_in_new$days_to_last_follow_up)
best.cluster=alist()
sig_autoen_var=alist()
for(k in 1:5){
  pval<-rep(0,35)
  for (i in 1:35){
    
    clin_in_new$cluster<-output[[k]]$fit[[i]]$clusters
    time=clin_in_new$overall_survival
    status=clin_in_new$deceased
    surv_diff <- survdiff(Surv(time, status)~cluster, data = clin_in_new)
    pval[i]<-surv_diff$pvalue
  }
  
  
  best.cluster[[k]]<-output[[k]]$fit[[which.min(pval)]]$clusters
  
  
  sig_autoen_var[[k]]<-output[[k]]$fit[[which.min(pval)]]$meanZ
  
}

#Harbasz Index
h<-c()
s<-c()
for(i in 1:5){
  harbasz<-function(x){
    return(calinhara(sig_autoen_var[[i]],x))
  }
  
  silhouette_score <- function(x){
    ss <- silhouette(x, dist(sig_autoen_var[[i]]))
    return(mean(ss[, 3]))
  }
  
  h[i]<-harbasz(best.cluster[[i]])
  s[i]<-silhouette_score(best.cluster[[i]])
  
}


#Harbarsz and silhoutte index plot
x<-seq(2,6)
clain_herb1<-h
sil1<-s

par(mar = c(5, 4,4 , 4) + 1.2)              # Additional space for second y-axis
plot(x, clain_herb1, col = "blue",type="b",ylab="Calinski-Harabasz")              # Create first plot

title("A")

par(new = TRUE)                             # Add new plot

plot(x, sil1, col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "Number of Clusters", ylab = "",type="b")
axis(side = 4, at = pretty(range(sil1)))      # Add second axis
mtext("Silhouette Index", side = 4, line = 4) 

legend("topright",legend=c("Calinski-Harabasz", "Silhouette Index"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"),cex=0.7,bty="n")



#Clustric metric visualization for sclae.lambda=c(0.10,0.65,0.92)

setwd("./Analysis_brca/icluster/cv_0.10_0.65_0.92/")

output=alist()
files=grep("cv.fit.icluster",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda) 
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)

setwd("C:/Users/pavel/Documents/")
clin_in_new<-read.table("./Analysis_brca/data/clin_brca_use.txt")
clin_in_new$deceased=clin_in_new$vital_status=="Dead"
clin_in_new$overall_survival=ifelse(clin_in_new$deceased,clin_in_new$days_to_death, clin_in_new$days_to_last_follow_up)
best.cluster=alist()
sig_autoen_var=alist()
for(k in 1:5){
  pval<-rep(0,35)
  for (i in 1:35){
    
    clin_in_new$cluster<-output[[k]]$fit[[i]]$clusters
    time=clin_in_new$overall_survival
    status=clin_in_new$deceased
    surv_diff <- survdiff(Surv(time, status)~cluster, data = clin_in_new)
    pval[i]<-surv_diff$pvalue
  }
  
  
  best.cluster[[k]]<-output[[k]]$fit[[which.min(pval)]]$clusters
  
  
  sig_autoen_var[[k]]<-output[[k]]$fit[[which.min(pval)]]$meanZ
  
}

#Harbasz Index
h<-c()
s<-c()
for(i in 1:5){
  harbasz<-function(x){
    return(calinhara(sig_autoen_var[[i]],x))
  }
  
  silhouette_score <- function(x){
    ss <- silhouette(x, dist(sig_autoen_var[[i]]))
    return(mean(ss[, 3]))
  }
  
  h[i]<-harbasz(best.cluster[[i]])
  s[i]<-silhouette_score(best.cluster[[i]])
  
}


#Harbarsz and silhoutte index plot
x<-seq(2,6)
clain_herb1<-h
sil1<-s

par(mar = c(5, 4,4 , 4) + 1.2)              # Additional space for second y-axis
plot(x, clain_herb1, col = "blue",type="b",ylab="Calinski-Harabasz")              # Create first plot

title("B")

par(new = TRUE)                             # Add new plot

plot(x, sil1, col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "Number of Clusters", ylab = "",type="b")
axis(side = 4, at = pretty(range(sil1)))      # Add second axis
mtext("Silhouette Index", side = 4, line = 4) 

legend("topright",legend=c("Calinski-Harabasz", "Silhouette Index"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"),cex=0.7,bty="n")

#Clustric metric visualization for sclae.lambda=c(1,0.2,1)

setwd("./Analysis_brca/icluster/")

output=alist()
files=grep("cv.ics_1_0.2_1",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda) 
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)

setwd("C:/Users/pavel/Documents/")
clin_in_new<-read.table("./Analysis_brca/data/clin_brca_use.txt")
clin_in_new$deceased=clin_in_new$vital_status=="Dead"
clin_in_new$overall_survival=ifelse(clin_in_new$deceased,clin_in_new$days_to_death, clin_in_new$days_to_last_follow_up)
best.cluster=alist()
sig_autoen_var=alist()
for(k in 1:5){
  pval<-rep(0,35)
  for (i in 1:35){
    
    clin_in_new$cluster<-output[[k]]$fit[[i]]$clusters
    time=clin_in_new$overall_survival
    status=clin_in_new$deceased
    surv_diff <- survdiff(Surv(time, status)~cluster, data = clin_in_new)
    pval[i]<-surv_diff$pvalue
  }
  
  
  best.cluster[[k]]<-output[[k]]$fit[[which.min(pval)]]$clusters
  
  
  sig_autoen_var[[k]]<-output[[k]]$fit[[which.min(pval)]]$meanZ
  
}

#Harbasz Index
h<-c()
s<-c()
for(i in 1:5){
  harbasz<-function(x){
    return(calinhara(sig_autoen_var[[i]],x))
  }
  
  silhouette_score <- function(x){
    ss <- silhouette(x, dist(sig_autoen_var[[i]]))
    return(mean(ss[, 3]))
  }
  
  h[i]<-harbasz(best.cluster[[i]])
  s[i]<-silhouette_score(best.cluster[[i]])
  
}


#Harbarsz and silhoutte index plot
x<-seq(2,6)
clain_herb1<-h
sil1<-s

par(mar = c(5, 4,4 , 4) + 1.2)              # Additional space for second y-axis
plot(x, clain_herb1, col = "blue",type="b",ylab="Calinski-Harabasz")              # Create first plot

title("C")

par(new = TRUE)                             # Add new plot

plot(x, sil1, col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "Number of Clusters", ylab = "",type="b")
axis(side = 4, at = pretty(range(sil1)))      # Add second axis
mtext("Silhouette Index", side = 4, line = 4) 

legend("topright",legend=c("Calinski-Harabasz", "Silhouette Index"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"),cex=0.7,bty="n")



#Clustric metric visualization for sclae.lambda=c(1,0.05,1)

setwd("./Analysis_brca/icluster/")

output=alist()
files=grep("cv.ics_1_0.05_1",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda) 
# nK = length(output)
# BIC = getBIC(output)
# devR = getDevR(output)

setwd("C:/Users/pavel/Documents/")
clin_in_new<-read.table("./Analysis_brca/data/clin_brca_use.txt")
clin_in_new$deceased=clin_in_new$vital_status=="Dead"
clin_in_new$overall_survival=ifelse(clin_in_new$deceased,clin_in_new$days_to_death, clin_in_new$days_to_last_follow_up)
best.cluster=alist()
sig_autoen_var=alist()
for(k in 1:5){
  pval<-rep(0,35)
  for (i in 1:35){
    
    clin_in_new$cluster<-output[[k]]$fit[[i]]$clusters
    time=clin_in_new$overall_survival
    status=clin_in_new$deceased
    surv_diff <- survdiff(Surv(time, status)~cluster, data = clin_in_new)
    pval[i]<-surv_diff$pvalue
  }
  
  
  best.cluster[[k]]<-output[[k]]$fit[[which.min(pval)]]$clusters
  
  
  sig_autoen_var[[k]]<-output[[k]]$fit[[which.min(pval)]]$meanZ
  
}

#Harbasz Index
h<-c()
s<-c()
for(i in 1:5){
  harbasz<-function(x){
    return(calinhara(sig_autoen_var[[i]],x))
  }
  
  silhouette_score <- function(x){
    ss <- silhouette(x, dist(sig_autoen_var[[i]]))
    return(mean(ss[, 3]))
  }
  
  h[i]<-harbasz(best.cluster[[i]])
  s[i]<-silhouette_score(best.cluster[[i]])
  
}


#Harbarsz and silhoutte index plot
x<-seq(2,6)
clain_herb1<-h
sil1<-s

par(mar = c(5, 4,4 , 4) + 1.2)              # Additional space for second y-axis
plot(x, clain_herb1, col = "blue",type="b",ylab="Calinski-Harabasz")              # Create first plot

title("D")

par(new = TRUE)                             # Add new plot

plot(x, sil1, col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "Number of Clusters", ylab = "",type="b")
axis(side = 4, at = pretty(range(sil1)))      # Add second axis
mtext("Silhouette Index", side = 4, line = 4) 

legend("topright",legend=c("Calinski-Harabasz", "Silhouette Index"),
       text.col=c("blue","red"),pch=c(15,15),col=c("blue","red"),cex=0.7,bty="n")



#Survival and Multidimensional Scalling plot

library(survival)
library(survminer)
library(ggplot2)


#Survival and Multidimensional Scalling plot for sclae.lambda=c(1,1,1)

clin_in1<-read.table("./Analysis_brca/icluster/clin_in_1_1_1_icls_brca.txt")


fit1<-survfit(Surv(overall_survival,deceased)~cls2_label,data=clin_in1)

print(fit1)
names(fit1$strata) <- gsub("cls2_label=", "", names(fit1$strata))
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
plot1<-ggsurvplot(fit1,data=clin_in1,legend = "right",pval=T,
                  title='Iclusterplus',
                  ggtheme=custom_theme(),legend.title="")

dist.matrix1 <- dist(clin_in1$cluster2_z)
MDS.2d1 <- cmdscale(dist.matrix1, k = 2)
options(repr.plot.width = 6, repr.plot.height=6)
#plot(MDS.2d[,1], MDS.2d[,2], col =integrated_clinical$lcluster5, pch = 16, main = "2D metric MDS Plot", xlab = "", ylab = "")
type1=as.factor(clin_in1$cls2_label)
p1<-as.data.frame((MDS.2d1)) %>% ggplot(aes(x=MDS.2d1[,1],y=MDS.2d1[,2],col=type1))+geom_point(size=3,fill=20)
plot2<-p1+theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Iclusterplus")+
  labs(color=NULL)

time=clin_in1$overall_survival
status=clin_in1$deceased
surv_diff <- survdiff(Surv(time, status)~cluster2, data = clin_in1)



#Survival and Multidimensional Scalling plot for sclae.lambda=c(0.10,0.65,0.92)


clin_in2<-read.table("./Analysis_brca/icluster/clin_in_cv_0.10_0.65_0.92_brca.txt")


fit2<-survfit(Surv(overall_survival,deceased)~cls2_label,data=clin_in2)

print(fit2)
names(fit2$strata) <- gsub("cls2_label=", "", names(fit2$strata))
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
plot3<-ggsurvplot(fit2,data=clin_in2,legend = "right",pval=T,
                  #  title='Autoencoder(100) with epoch=50', 
                  ggtheme=custom_theme(),legend.title="")


dist.matrix2 <- dist(clin_in2$cluster2_z)
MDS.2d2 <- cmdscale(dist.matrix2, k = 2)
options(repr.plot.width = 6, repr.plot.height=6)
#plot(MDS.2d[,1], MDS.2d[,2], col =integrated_clinical$lcluster5, pch = 16, main = "2D metric MDS Plot", xlab = "", ylab = "")
type2=as.factor(clin_in2$cls2_label)
p2<-as.data.frame((MDS.2d2)) %>% ggplot(aes(x=MDS.2d2[,1],y=MDS.2d2[,2],col=type2))+geom_point(size=3,fill=20)
plot4<-p2+theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                plot.title = element_text(hjust = 0.5))+
  #ggtitle("Autoencoder(100) with epoch=50")+
  labs(color=NULL)

time=clin_in2$overall_survival
status=clin_in2$deceased
surv_diff <- survdiff(Surv(time, status)~cluster2, data = clin_in2)


#Survival and Multidimensional Scalling plot for sclae.lambda=c(1,0.2,1)


clin_in3<-read.table("./Analysis_brca/icluster/clin_in_1_0.2_1_icls_brca.txt")


fit3<-survfit(Surv(overall_survival,deceased)~cls2_label,data=clin_in3)

print(fit3)
names(fit3$strata) <- gsub("cls2_label=", "", names(fit3$strata))
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
plot5<-ggsurvplot(fit3,data=clin_in3,legend = "right",pval=T,
                  # title='Autoencoder(150) with epoch=50', 
                  ggtheme=custom_theme(),legend.title="")


dist.matrix3 <- dist(clin_in3$cluster2_z)
MDS.2d3 <- cmdscale(dist.matrix3, k = 2)
options(repr.plot.width = 6, repr.plot.height=6)
#plot(MDS.2d[,1], MDS.2d[,2], col =integrated_clinical$lcluster5, pch = 16, main = "2D metric MDS Plot", xlab = "", ylab = "")
type3=as.factor(clin_in3$cls2_label)
p3<-as.data.frame((MDS.2d3)) %>% ggplot(aes(x=MDS.2d3[,1],y=MDS.2d3[,2],col=type3))+geom_point(size=3,fill=20)
plot6<-p3+theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                plot.title = element_text(hjust = 0.5))+
  # ggtitle("Autoencoder(150) with epoch=50")+
  labs(color=NULL)

time=clin_in3$overall_survival
status=clin_in3$deceased
surv_diff <- survdiff(Surv(time, status)~cluster2, data = clin_in3)


#Survival and Multidimensional Scalling plot for sclae.lambda=c(1,0.05,1)


clin_in4<-read.table("./Analysis_brca/icluster/clin_in_1_0.05_1_icls_brca.txt")


fit4<-survfit(Surv(overall_survival,deceased)~cls2_label,data=clin_in4)

print(fit4)
names(fit4$strata) <- gsub("cls2_label=", "", names(fit4$strata))
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

plot7<-ggsurvplot(fit4,data=clin_in4,pval=T,legend = "right",
                  #   title='Autoencoder(200) with epoch=50', 
                  ggtheme=custom_theme(),legend.title="")


dist.matrix4 <- dist(clin_in4$cluster2_z)
MDS.2d4 <- cmdscale(dist.matrix4, k = 2)
options(repr.plot.width = 6, repr.plot.height=6)
#plot(MDS.2d[,1], MDS.2d[,2], col =integrated_clinical$lcluster5, pch = 16, main = "2D metric MDS Plot", xlab = "", ylab = "")
type4=as.factor(clin_in4$cls2_label)
p4<-as.data.frame((MDS.2d4)) %>% ggplot(aes(x=MDS.2d4[,1],y=MDS.2d4[,2],col=type4))+geom_point(size=3,fill=20)
plot8<-p4+theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                plot.title = element_text(hjust = 0.5))+
  # ggtitle("Autoencoder(200) with epoch=50")+
  labs(color=NULL)

time=clin_in4$overall_survival
status=clin_in4$deceased
surv_diff <- survdiff(Surv(time, status)~cluster2, data = clin_in4)

# splot1<-ggarrange(plot1$plot,plot2,plot3$plot,plot4,plot5$plot,plot6,plot7$plot,plot8,
#                   labels=c("A","B","C","D","E","F","G","H"),ncol=2,nrow=4)


splot1<-ggarrange(plot1$plot,plot2,plot3$plot,plot4,labels=c("A","B","C","D"),ncol=2,nrow=2)

splot2<-ggarrange(plot5$plot,plot6,plot7$plot,plot8,labels=c("E","F","G","H"),ncol=2,nrow=2)

