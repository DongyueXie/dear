#'@param Y sample by feature matrix
#'@param X sample by covriates design matrix
#'@param n.sv number of surrogate variables
#'@param B number of iterations in sva

slash = function(Y,X,n.sv=NULL,B=5){
  #apply sva
  sva_sva = irwsva.build(t(Y),X[,-1,drop=FALSE],X[,1,drop=FALSE],n.sv,B)
  X.lm = cbind(X[,1]-rowSums(X[,-1,drop=FALSE]),X[,-1])
  X.lm = cbind(X.lm,sva_sva$sv)
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
  if(length(strong)>1){
    U.pca = cov_pca(mash_data.L,2,subset=strong)
    U.ed = cov_ed(mash_data.L, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data.L, c(U.c,U.ed), algorithm.version = 'Rcpp')
  }else{
    m = mash(mash_data.L, c(U.c), algorithm.version = 'Rcpp')
  }

  return(m)
}



#'@param Y sample by feature matrix
#'@param X sample by covriates design matrix

lash = function(Y,X){
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
  if(length(strong)>1){
    U.pca = cov_pca(mash_data.L,2,subset=strong)
    U.ed = cov_ed(mash_data.L, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data.L, c(U.c,U.ed), algorithm.version = 'Rcpp')
  }else{
    m = mash(mash_data.L, c(U.c), algorithm.version = 'Rcpp')
  }

  return(m)
}


#'@title create signals for normalized data.
#'@param Y cell by gene matrix, normalized.
#'@param x a vector treatment assignments.
#'@param null_gene_prop proporiton of null genes
#'@param null_trt_prop propotion of null treatment
#'@param eff_sd sd of effect from normal distribution

create_signal = function(Y,x,null_gene_prop=0.9,null_trt_prop=0.8,eff_sd = 0.8){
  N = nrow(Y)
  G = ncol(Y)
  Ntrt = length(unique(x))
  non_null_gene_idx = sample(1:G,round(G*(1-null_gene_prop)),replace = FALSE)

  non_null_trt_idx = matrix(nrow=round(Ntrt*(1-null_trt_prop)),ncol=G)
  for(g in non_null_gene_idx){
    non_null_trt = sample(1:Ntrt,round(Ntrt*(1-null_trt_prop)),replace = FALSE)
    non_null_trt_idx[,g] = non_null_trt
    non_null_trt_effect = c()
    for(i in non_null_trt[-1]){
      idx = which(x==i)
      effect = rnorm(1,0,eff_sd)
      Y[idx,g] = Y[idx,g] + effect
      non_null_trt_effect = c(non_null_trt_effect,effect)
    }
    idx = which(x==non_null_trt[1])
    Y[idx,g] = Y[idx,g] - sum(non_null_trt_effect)
  }
  return(list(Y=Y,non_null_gene_idx = non_null_gene_idx, non_null_trt_idx=non_null_trt_idx,group_idx = x))
}


