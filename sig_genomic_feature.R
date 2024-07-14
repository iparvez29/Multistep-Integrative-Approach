# Maintained and edited by Imran Parvez, email id: iparvez@augusta.edu
# This code assumed that the user has downloaded the three types of 
# genomic data Rna expression, CNV and DNA Methylation from TCGA database.

# The all three types of data has cleaned by deleting all the features 
# having at least one missing value or 80% of samples having zero values.
# The datasets has to matched across all three types of data.

# The datasets were organized as subjects are columns
# rows are features corresponding to each subjects.

# This codes are used to find the significant features in univariate cox-PH model.

# Finding the significant Rna expression features which significantly contribute 
# to the Overall Survival of Breast Cancer Patients

# This codes were run in Augusta University High Performance Computing (HPC) Cluster. 

# Load required packages
library(survival)
library(parallel)

# The'clin_brca' contains survival data for breast cancer patients which are
# loaded from the specified path.

clin_brca<-read.table("/scratch/iparvez/analysis_brca/data/clin_brca_use.txt")

# The exp_nonazro_brca contains the Rna Expression for breast cancer patients.
exp_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/exp_nonazro_brca.txt")
#The data has transposed so that the rows are subjects and the columns are features.
exp_nonazro_brca<-t(exp_nonazro_brca)
exp_nonazro_brca<-as.data.frame(exp_nonazro_brca)


rna_data<-cbind(clin_brca$deceased,clin_brca$overall_survival,exp_nonazro_brca)
data<-rna_data
cox<-function(v) {
  
  cx<-summary(coxph(Surv(clin_brca$overall_survival,clin_brca$deceased)~v,data=data))
  return(cx$sctest["pvalue"])
}

cl<-makeCluster(20)
clusterExport(cl,c("data","cox","clin_brca"))
clusterEvalQ(cl, {library(survival)})


exp_pvalue<-parApply(cl,data[,-c(1:2)],2,FUN=cox)
# The Rna expression features selected for p-value<001
exp_sig<-as.matrix(exp_nonazro_brca[,exp_pvalue<0.01])
stopCluster(cl)
#The significant Rna expression features are saved in to a text file "exp_sig_icluster_p0.01.txt"
# in the specified path
write.table(exp_sig,"/scratch/iparvez/analysis_brca/data/exp_sig_icluster_p0.01.txt")


# Finding the significant Copy Number Variation(CNV) features which significantly
# contribute to the Overall Survival of Breast Cancer Patients

#The CNV data has read from a specified path
copy_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/copy_nonazro_brca.txt")

#The data has transposed so that the rows are subjects and the columns are features.
copy_nonazro_brca<-t(copy_nonazro_brca)
copy_nonazro_brca<-as.data.frame(copy_nonazro_brca)

cp_data<-cbind(clin_brca$deceased,clin_brca$overall_survival,copy_nonazro_brca)
data<-cp_data

cl<-makeCluster(20)
clusterExport(cl,c("data","cox","clin_brca"))
clusterEvalQ(cl, {library(survival)})

pvalue_cp<-parApply(cl,cp_data[,-c(1:2)],2,FUN=cox)

# The CNV features selected for p-value<001
sig_copy<-as.matrix(copy_nonazro_brca[,pvalue_cp<0.01])
stopCluster(cl)


#The significant CNV features are saved in to a text file "sig_copy_icluster_p0.01.txt"
# in the specified path

write.table(sig_copy,"/scratch/iparvez/analysis_brca/data/sig_copy_icluster_p0.01.txt")



# Finding the significant DNA Methylation features which significantly 
# contribute to the Overall Survival of Breast Cancer Patients

#The methylatin data has read from a specified path
methylm_nonazro_brca<-read.table("/scratch/iparvez/analysis_brca/data/methylm_nonazro_brca.txt")
#The data has transposed so that the rows are subjects and the columns are features.
methylm_nonazro_brca<-t(methylm_nonazro_brca)

methylm_nonazro_brca<-as.data.frame(methylm_nonazro_brca)

methyl_data<-cbind(clin_brca$deceased,clin_brca$overall_survival,methylm_nonazro_brca)
data<-methyl_data

cl<-makeCluster(20)
clusterExport(cl,c("data","cox","clin_brca"))
clusterEvalQ(cl, {library(survival)})

methyl_pvalue<-parApply(cl,data[,-c(1:2)],2,cox)

# The methylation features selected for p-value<001
methyl_sig<-as.matrix(methylm_nonazro_brca[,methyl_pvalue<0.01])

stopCluster(cl)

#The significant methylation features are saved in to a text file 
# "sig_copy_icluster_p0.01.txt" in the specified path
write.table(methyl_sig,"/scratch/iparvez/analysis_brca/data/methyl_sig_icluster_p0.01.txt")



