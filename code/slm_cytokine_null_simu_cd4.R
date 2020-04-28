source('code/slm.R')
library(mashr)
library(Matrix)
library(sva)
library(RSpectra)
load("data/cytokine_cd4.RData")
#remove genes that appears in less than 50*5 cells
ngroup = 50
rm.idx = which(colSums(data.cd4!=0)<(ngroup*5))
#rm.idx = which(colSums(data.cd4)==0)
data.cd4 = data.cd4[,-rm.idx]
data.cd4 = as.matrix(data.cd4)


ncell = nrow(data.cd4)
ngene = ncol(data.cd4)

n_simu = 10
out_slm = list()
out_lm = list()

set.seed(12345)
for(i in 1:n_simu){
  print(paste("running simulation: ", i))
  group_idx = sample(1:ngroup,ncell,replace = TRUE)
  X = model.matrix(~as.factor(group_idx))
  out_slm[[i]] = slash(data.cd4,X,n.sv=20,B=2)
  save(out_slm,file='output/slm_cytokine_null_simu_cd4.RData')
  out_lm[[i]] = lash(data.cd4,X)
  save(out_lm,file='output/lm_cytokine_null_simu_cd4.RData')
}
