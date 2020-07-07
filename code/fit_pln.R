
# fit pln to bulk_cytokine

load('D:/dear/data/bulk_cytokin_SI.RData')
dim(mat)
samples = colnames(mat)
samples = matrix(unlist(strsplit(samples,split = '_')),ncol=3,byrow = TRUE)

# filter out genes that 1. have at least one zero count in any sample; 2. the union of top 1% expressed gene in each condition
library(Matrix)
mat = mat[which(rowSums(mat==0)<1),]

rm.idx = c()
for(i in 1:ncol(mat)){

  rm.idx.i1 = which(mat[,i]>=quantile(mat[,i],0.95))
  rm.idx.i2 = which(mat[,i]<=quantile(mat[,i],0.05))
  rm.idx = c(rm.idx,rm.idx.i1,rm.idx.i2)

}

rm.idx = unique(rm.idx)

mat = mat[-rm.idx,]

library(PLNmodels)
colnames(mat) = 1:ncol(mat)
datax = prepare_data(counts = t(mat[,1:30]),covariates = data.frame(group=samples[1:30,2]),offset = 'RLE')
myPLN_full <- PLN(Abundance ~ 1+offset(log(Offset))+group, datax,control = list(covariance = "full", trace = 0))
save(myPLN_full,file='output/PLNfit_bulk_cytokin_SI_1_30.RData')
rm(myPLN_full)

datax = prepare_data(counts = t(mat[,31:90]),covariates = data.frame(group=samples[31:90,2]),offset = 'RLE')
myPLN_full <- PLN(Abundance ~ 1+offset(log(Offset))+group, datax,control = list(covariance = "full", trace = 0))
save(myPLN_full,file='output/PLNfit_bulk_cytokin_SI_31_90.RData')
