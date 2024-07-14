
# This code uses the significant features of RNA expression, CNV and DNA methylation
# selected from the R-code shared named as "sig_genomic_feature.R".

# The datasets were organized as subjects are rows
# columns are features corresponding to each subjects.

# Loaded required packages
library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)

# The directory are set to the specified path
setwd("C:/Users/pavel/Documents/Analysis_brca/icluster/")

# The'clin_brca' contains survival data for breast cancer patients which are
# loaded from the specified path.
clin_in<-read.table("./Analysis_brca/data/clin_brca_use.txt")

# This function is defined for mean absolute deviation (MAD+scale) normalization.

mad_nrmlz<-function(x){
  
  mad<-median(as.numeric(abs(x-median(as.numeric(x)))))
  return((x-median(as.numeric(x)))/(mad+1e-4))
}


# This function is defined for Z-score normalization.

mean_nrmlz<-function(x){
  
  sd<-sqrt(var(x))
  return((x-mean(as.numeric(x)))/sd)
}


# This function is defined for robust scale normalization.

robust_nrmlz<-function(x){
  x<-as.numeric(x)
  m<-mean(x[x>quantile(x,0.25) & x<quantile(x,0.75)])
  s<-sqrt(var(x[x>quantile(x,0.25) & x<quantile(x,0.75)]))
  return((x-m)/s)
  
}


# The all three types of significant genomic features has downloaded from 
# the specified path which were found from the shared R-code "sig_genomic_feature.R".

# The exp_sig is the dataset contains significant Rna expression features.
exp_sig<-read.table("./Analysis_brca/icluster/exp_sig_icluster_p0.01.txt")

# The sig_copy is the dataset contains significant CNV features.
sig_copy<-read.table("./Analysis_brca/icluster/sig_copy_icluster_p0.01.txt")

# The methyl_sig is the data sets contains significant methylation features
methyl_sig<-read.table("./Analysis_brca/icluster/methyl_sig_icluster_p0.01.txt")

# All three types of significant genomic features are normalized 
# by the corresponding function
exp_sig<-apply(exp_sig,2,mad_nrmlz)
sig_copy<-apply(sig_copy,2,mean_nrmlz)
methyl_sig<-apply(methyl_sig,2,robust_nrmlz)

# The integrative clustering has performed by using iCluster+ method.
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=10,dt1=exp_sig,dt2=sig_copy,dt3=methyl_sig,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,n.lambda=35, scale.lambda=c(1,1,1),maxiter=20)
  save(cv.fit, file=paste("./Analysis_brca/icluster/cv_ics_1_1_1_brca.k",k,".Rdata",sep=""))
}

# The results found from the iCluster+ has listed in to the list variable output
output=alist()
files=grep("cv.ics_1_1_1",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda) 


setwd("C:/Users/pavel/Documents/")

clin_in_new<-read.table("./Analysis_brca/data/clin_brca_use.txt")
clin_in_new$deceased=clin_in_new$vital_status=="Dead"
clin_in_new$overall_survival=ifelse(clin_in_new$deceased,clin_in_new$days_to_death, clin_in_new$days_to_last_follow_up)

# The pvalue of survival difference has computed by using survival difference test
# for all the clusters found from iCluster+. 
# Also the best cluster has selected for each k (ranges from 1 to 5) values
# based on smallest p-value.

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


