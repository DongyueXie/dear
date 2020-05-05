source('code/slm.R')
library(mashr)
library(Matrix)
library(sva)
library(RSpectra)
load("data/cytokine_cd4.RData")
#remove genes that appears in less than 50*5 cells
ngroup = 50
rm.idx = which(colSums(data.cd4!=0)<(ngroup*20))
#rm.idx = which(colSums(data.cd4)==0)
data.cd4 = data.cd4[,-rm.idx]
Y = as.matrix(data.cd4)

nonzero.elements = c(Y)[which(c(Y)!=0)]




ncell = nrow(Y)
ngene = ncol(Y)

set.seed(12345)
group_idx = sample(1:ngroup,ncell,replace = TRUE)
X = model.matrix(~as.factor(group_idx))



#linear model
X.lm = cbind(X[,1]-rowSums(X[,-1,drop=FALSE]),X[,-1])
cov_of_interest = 1:ncol(X)
lmout = limma::lmFit(object = t(Y), design = X.lm)
eout  = limma::eBayes(lmout)
out = list()
out$betahat = lmout$coefficients[,cov_of_interest,drop=FALSE]
out$sebetahat = lmout$stdev.unscaled[,cov_of_interest,drop=FALSE] * sqrt(eout$s2.post)


mash_data = mash_set_data(out$betahat,out$sebetahat)
mash_data.L = mash_update_data(mash_data,ref='mean')
# canonical cov of mash
U.c = cov_canonical(mash_data.L)
# data-driven cov of mash
m.1by1 = mash_1by1(mash_data.L)
strong = get_significant_results(m.1by1)
U.pca = cov_pca(mash_data.L,2,subset=strong)
U.ed = cov_ed(mash_data.L, U.pca, subset=strong)
# fit mash model
m = mash(mash_data.L, c(U.c,U.ed), algorithm.version = 'R')


#yusha

extract.data <- function(data,annot,trt){

  Bhat <- matrix(NA,nrow=nrow(data),ncol=length(trt))
  Shat <- matrix(NA,nrow=nrow(data),ncol=length(trt))
  rownames(Bhat) <- rownames(data)
  colnames(Bhat) <- trt
  rownames(Shat) <- rownames(data)
  colnames(Shat) <- trt

  overall.mean <- Matrix::rowMeans(data)

  for(i in 1:length(trt)){
    data.tmp <- data[, annot==trt[i]]
    trt.mean <- Matrix::rowMeans(data.tmp)
    Bhat[,i] <- trt.mean
    Shat[,i] <- sqrt((apply(data.tmp, 1, var)+(trt.mean-overall.mean)^2)/ncol(data.tmp))
  }

  return(list(Bhat=Bhat, Shat=Shat))
}

out.s = extract.data(t(Y),group_idx,1:50)


mash_data.s = mash_set_data(out.s$Bhat,out.s$Shat)
mash_data.L.s = mash_update_data(mash_data.s,ref='mean')
# canonical cov of mash
U.c.s = cov_canonical(mash_data.L.s)
# data-driven cov of mash
m.1by1.s = mash_1by1(mash_data.L.s)
strong.s = get_significant_results(m.1by1.s)

U.pca.s = cov_pca(mash_data.L.s,2,subset=strong.s)
U.ed.s = cov_ed(mash_data.L.s, U.pca.s, subset=strong.s)
# fit mash model
m.s = mash(mash_data.L.s, c(U.c.s,U.ed.s), algorithm.version = 'R')

### compare

compare_lm_sample = list(Y=Matrix(Y,sparse=TRUE),rm.idx=rm.idx,group_idx=group_idx,out_lm = out, mash_lm = m, out_sam = out.s,mash_sam=m.s)

save(compare_lm_sample,file='output/compare_lm_sample.RData')


as.numeric(out$betahat[279,])
as.numeric(out.s$Bhat[279,])

boxplot(cbind(as.numeric(out$betahat[279,]),as.numeric(out.s$Bhat[279,])))

as.numeric(out$sebetahat[279,])
as.numeric(out.s$Shat[279,])

boxplot(as.numeric(out$sebetahat[279,]),as.numeric(out.s$Shat[279,]))
