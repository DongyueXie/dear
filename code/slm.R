#'@param Y sample by feature matrix
#'@param X sample by covriates design matrix
#'@param n.sv number of surrogate variables
#'@param B number of iterations in sva

slash = function(Y,X,n.sv=NULL,B=5){
  #apply sva
  sva_sva = irwsva.build(t(Y),X[,-1,drop=FALSE],X[,1,drop=FALSE],n.sv,B)
  X.lm = cbind(X[,-1],X[,1]-rowSums(X[,-1,drop=FALSE]))
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
  X.lm = cbind(X[,-1],X[,1]-rowSums(X[,-1,drop=FALSE]))
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

