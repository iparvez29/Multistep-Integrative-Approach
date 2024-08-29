# This code assumed that the user has downloaded the three types of 
# genomic data Rna expression, CNV and DNA Methylation from TCGA database.

# The all three types of data has cleaned by deleting all the features 
# having at least one missing value or 80% of samples having zero values.
# The datasets has to matched across all three types of data.

# The datasets were organized as subjects are columns
# rows are features corresponding to each subjects.


# Loaded the required packages
library(dplyr)
library(caret)
library(sparsegl)
library(combinat)
library(doParallel)


# Rna expression data has read from the specified path
exp_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/exp_nonazro_brca.txt")
# The Rna expression data has been transposed so that each row represents 
# a subject, and each column corresponds to a feature for those subjects.
exp_nonazro_brca<-t(exp_nonazro_brca)
# CNV data has read from the specified path
copy_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/copy_nonazro_brca.txt")
# The CNV data has been transposed so that each row represents a subject,
# and each column corresponds to a feature for those subjects.
copy_nonazro_brca<-t(copy_nonazro_brca)
# DNA methylation data has read from the specified path
methylm_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/methylm_nonazro_brca.txt")
# The methylation data has been transposed so that each row represents a subject,
# and each column corresponds to a feature for those subjects.

methylm_nonazro_brca<-t(methylm_nonazro_brca)

# The'clin_in' contains survival data for breast cancer patients which are
# loaded from the specified path.
clin_in<-read.table("/scratch/iparvez/analysis_brca/icluster/clin_in_1_0.05_1_icls_brca.txt")

# The wilcoxon function is defined to perform wilcoxon rank sum test between
# two survival groups (high survival and low survival) and identify the 
# significantly differentially expressed features between these two groups.
wilcoxon<-function(x){
  
  w<-wilcox.test(x~clin_in$cluster2,
                 exact = FALSE)
  return(w$p.value)
}


# The wilcoxon rank sum test has been performed on Rna expression
exp_pvalue<-apply(exp_nonazro_brca,2,wilcoxon)
p_adjst_exp<-p.adjust(exp_pvalue,method="fdr")
exp_sigvar<-as.matrix(exp_nonazro_brca[,(which(p_adjst_exp<0.05))])

# The wilcoxon rank sum test has been performed on CNV data
copy_pvalue<-apply(copy_nonazro_brca,2,wilcoxon)
p_adjst_cp<-p.adjust(copy_pvalue,method="fdr")
copy_sigvar<-as.matrix(copy_nonazro_brca[,(which(p_adjst_cp<0.05))])
# The wilcoxon rank sum test has been performed on Dna metyhlation data
methyl_pvalue<-apply(methylm_nonazro_brca,2,wilcoxon)
p_adjst_methyl<-p.adjust(methyl_pvalue,method="fdr")
methyl_sigvar<-as.matrix(methylm_nonazro_brca[,(which(p_adjst_methyl<0.05))])


group<-c(rep(1,dim(exp_sigvar)[2]),rep(2,dim(copy_sigvar)[2]),rep(3,dim(methyl_sigvar)[2]))


cancer_type<-substr(clin_in$project_id,start=6,stop=9)

#Th fld function is defined to split the dataset into train dataset and test dataset. 
fld<-function(x) {
  set.seed(1234)
  idx <- createFolds(as.factor(cancer_type), k=5)
  trn<-idx[-c(x[1],x[2])]
  train<-unlist(trn)
  test<-unlist(idx[c(x[1],x[2])])
  return(list(train=train,test=test))
}
cmbn<-as.matrix(combn(5,2))

trn_tst<-apply(cmbn,2,fld)

#These three functions are defined to normailize the three types of genomic data.

mad_nrmlz<-function(x){
  
  mad<-median(as.numeric(abs(x-median(as.numeric(x)))))
  return((x-median(as.numeric(x)))/mad)
}


mean_nrmlz<-function(x){
  
  sd<-sqrt(var(x))
  return((x-mean(as.numeric(x)))/sd)
}


robust_nrmlz<-function(x){
  x<-as.numeric(x)
  m<-mean(x[x>quantile(x,0.25) & x<quantile(x,0.75)])
  s<-sqrt(var(x[x>quantile(x,0.25) & x<quantile(x,0.75)]))
  return((x-m)/s)
  
}



y0<-clin_in$cluster2



exp_nrmlz_trn<-apply(exp_sigvar,2,mad_nrmlz)
copy_nrmlz_trn<-apply(copy_sigvar,2,mean_nrmlz)
methyl_nrmlz_trn<-apply(methyl_sigvar,2,robust_nrmlz)
x_trn<-cbind(exp_nrmlz_trn,copy_nrmlz_trn, methyl_nrmlz_trn)


fit_trial <- sparsegl(x_trn,y0,group, family = "binomial")
lmbda_given<-fit_trial$lambda

gc()

start_time<-Sys.time()
# numCores <- detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
# The ten-fold cross validation
mserr<-c()
clusterExport(cl,c("trn_tst","y0","mad_nrmlz","mean_nrmlz","robust_nrmlz","exp_sigvar","copy_sigvar","methyl_sigvar","mserr","group","lmbda_given"),envir=environment())
output<-foreach(i=1:10,.combine="cbind",.packages=c("sparsegl","caret")) %dopar% {
  
  exp_nrmlz_train<-apply(exp_sigvar[trn_tst[[i]]$train,],2,mad_nrmlz)
  copy_nrmlz_train<-apply(copy_sigvar[trn_tst[[i]]$train,],2,mean_nrmlz)
  methyl_nrmlz_train<-apply(methyl_sigvar[trn_tst[[i]]$train,],2,robust_nrmlz)
  x_train<-cbind(exp_nrmlz_train,copy_nrmlz_train, methyl_nrmlz_train)
  
  exp_nrmlz_test<-apply(exp_sigvar[trn_tst[[i]]$test,],2,mad_nrmlz)
  copy_nrmlz_test<-apply(copy_sigvar[trn_tst[[i]]$test,],2,mean_nrmlz)
  methyl_nrmlz_test<-apply(methyl_sigvar[trn_tst[[i]]$test,],2,robust_nrmlz)
  x_test<-cbind( exp_nrmlz_test,copy_nrmlz_test,  methyl_nrmlz_test)
  
  fit_logit <- sparsegl(x_train, y0[trn_tst[[i]]$train], group, family = "binomial",lambda = lmbda_given)
  prediction<-predict(fit_logit,new=x_test,s=fit_logit$lambda,type="class")
  prd.matrix<-matrix(prediction,length(y0[trn_tst[[i]]$test]),length(fit_logit$lambda))
  
  error<-function(py){
    py<-as.factor(py)
    
    
    tbl<-confusionMatrix(py,as.factor(y0[trn_tst[[i]]$test]))$table
    
    mce<-(tbl[1,2]+tbl[2,1])/length(y0[trn_tst[[i]]$test])
    
    
    return(mce)
  }
  
  
  
  mserr<-apply(prd.matrix,2,error)
  return(cbind(fit_logit$lambda,mserr))
  gc() 
}




stopCluster(cl)

write.table(output,"/scratch/iparvez/analysis_brca/icluster/out_icls_1_0.05_1_brca.txt")

