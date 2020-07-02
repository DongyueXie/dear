
library(seqgendiff)



lash = function(vout,mash_ref = "mean"){
  lmout = limma::lmFit(vout)
  out = list()
  out$betahat = lmout$coefficients
  out$sebetahat = lmout$stdev.unscaled * sqrt(eout$s2.post)
  mash_data = mash_set_data(out$betahat,out$sebetahat)
  mash_data.L = mash_update_data(mash_data,ref=mash_ref)
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

# @description Compute approximately unbiased variance estimates for
#   the estimators for logit(p) when n is small.
# @param n number of trials
# @param s number of successes
# @param f number of failures
v3 = function (n, s, f)
  return((n + 1)/n * (1/(s + 1) + 1/(f + 1)))

# @description Compute approximately unbiased variance estimates for
#   the estimators for logit(p) when n is small.
# @param n number of trials
# @param s number of successes
# @param f number of failures
vs = function (n, s, f) {
  vv = v3(n, s, f)
  return(vv * (1 - 2/n + vv/2))
}

# @description Compute approximately unbiased variance estimates for
#   the estimators for logit(p) when n is small.
# @param n number of trials
# @param s number of successes
# @param f number of failures
vss = function (n, s, f) {
  vv = v3(n, s, f)
  return(vs(n, s, f) - 1/2 * vv^2 * (vv - 4/n))
}

# X: cell by genes

pseudo_bulk = function(X,groups){
  n_group = length(unique(groups))
  n_gene = ncol(X)
  group_name = unique(groups)
  Xa = matrix(nrow=n_group,ncol=n_gene)
  for(i in 1:n_group){
    group_idx = which(groups==group_name[i])
    Xa[i,] = colSums(X[group_idx,])
  }
  rownames(Xa) = group_name
  Xa
}


################# ############################# ###############
################# mash on single cell real data ###############

log_ratio_approx = function(x,n){
  if(x==0){
    a=log((x+1/2)/(n-x+1/2)) - 1/2
  }else if(x==n){
    a=log((x+1/2)/(n-x+1/2)) + 1/2
  }else{
    a=log((x+1/2)/(n-x+1/2))
  }
  v = vss(n,x,n-x)
  return(list(a=a,v=v))
}

log_ratio_approx_cov = function(a1,a2,n1,n2,x0){
  x0*(1+exp(a1))*(1+exp(a2))/(n1*n2)
}

log_ratio_approx_vector = function(x,ref_x){
  n = length(x)
  Sigma = matrix(nrow=n,ncol=n)
  alpha = c()
  for(i in 1:n){
    out = log_ratio_approx(x[i],ref_x+x[i])
    alpha[i] = out$a
    Sigma[i,i] = out$v
  }

  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Sigma[i,j] = log_ratio_approx_cov(alpha[i],alpha[j],ref_x+x[i],ref_x+x[j],ref_x)
    }
  }
  Sigma = Sigma + t(Sigma) - diag(Sigma)

  return(list(alpha=alpha,Sigma=Sigma))
}


#'@title
#'

# mat: cell by genes

log_ratio_est = function(mat,groups,ref=1){
  mat_bulk = pseudo_bulk(mat,groups)
  nc = length(unique(groups))
  ng = ncol(mat)
  log_ratio = matrix(nrow=ng,ncol=nc)
  log_ratio_cov  = array(dim=c(nc,nc,ng))
  ref_group = which(unique(groups)==ref)
  for(i in 1:ng){
    out = log_ratio_approx_vector(mat_bulk[-ref_group,i],mat_bulk[ref_group,i])
    log_ratio[i,] = out$alpha
    log_ratio_cov[,,i] = out$Sigma
  }

  return(list(log_ratio = log_ratio,log_ratio_cov=log_ratio_cov))

}







# read data
library(Matrix)
samples = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_samples.rds')
datax = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_counts.rds')
cell = "CD4_T_cells"

cell_idx = which(samples$cell_type == cell)
data.cd4 = datax[cell_idx,]
data.cd4 = as.matrix(data.cd4)
# remove genes
ngroup = 30
rm.idx = which(colSums(data.cd4!=0)<(ngroup*5))
#rm.idx = which(colSums(data.cd4)==0)
data.cd4 = data.cd4[,-rm.idx]


# how many non-null genes
non_null_gene_prop = 0.1
non_null_trt_prop = 0.2

signal_sd = 0.8

out_null_voom = list()
out_null_agg = list()
#out_null_glm = list()

out_nonnull_voom = list()
out_nonnull_agg = list()
#out_nonnull_glm = list()

n_simu = 10
for(simu in 1:n_simu){
  # subset your data so that results do not depend on the quirks of a few individuals/genes.
  submat = t(select_counts(mat = t(data.cd4), nsamp = ngroup*300, ngene = 2000))

  rm.idx = which(colSums(submat!=0)<(ngroup*5))
  submat = submat[,-rm.idx]
  # generate design matrix
  ncell = nrow(submat)
  group_idx = sample(1:ngroup,ncell,replace = TRUE)
  X = model.matrix(~0+as.factor(group_idx))

  # Null simulation

  ## aggregate data approach
  out = log_ratio_est(submat,group_idx,ref = 1)
  mash_data = mash_set_data(out$log_ratio,V = out$log_ratio_cov)
  #mash_data.L = mash_update_data(mash_data,ref=mash_ref)
  # canonical cov of mash
  U.c = cov_canonical(mash_data)
  # data-driven cov of mash
  m.1by1 = mash_1by1(mash_data)
  strong = get_significant_results(m.1by1)
  if(length(strong)>1){
    U.pca = cov_pca(mash_data,2,subset=strong)
    U.ed = cov_ed(mash_data, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data, c(U.c,U.ed), algorithm.version = 'Rcpp')
  }else{
    m = mash(mash_data, c(U.c), algorithm.version = 'Rcpp')
  }
  out_null_agg[[simu]] = m
  save(out_null_agg,file = '/scratch/midway2/dyxie/sc-cytokine/output/cd4_null_pseudobulk.RData')


  ## directly use linear model
  vout = limma::voom(counts = t(submat), design = X)
  fit_null_voom = lash(vout,ref="1")
  out_null_voom[[simu]] = fit_null_voom

  save(out_null_voom,file = '/scratch/midway2/dyxie/sc-cytokine/output/cd4_null_voom.RData')


  # ## directly use generalized linear model
  # dds = DESeqDataSetFromMatrix(countData = t(submat),
  #                               colData = group_idx,
  #                               design = ~ 0+group_idx)
  # fit_null_glm =


  # Thin data

  designmat = X[,-1]


  coefmat = matrix(stats::rnorm(ncol(designmat) * ncol(submat),0,signal_sd),
                    ncol = ncol(designmat),
                    nrow = ncol(submat))

  non_null_matrix = matrix(0,nrow = ncol(submat),ncol = ncol(designmat))
  for(j in sample(1:ncol(submat),round(ncol(submat)*non_null_gene_prop))){
    non_null_matrix[j,sample(1:ncol(designmat),round(ncol(designmat)*non_null_trt_prop))] = 1
  }

  coefmat = coefmat*non_null_matrix

  thout = thin_diff(mat = t(submat), design_fixed = designmat,  coef_fixed   = coefmat)

  ##
  ## aggregate data approach
  out = log_ratio_est(t(thout$mat),group_idx,ref = 1)
  mash_data = mash_set_data(out$log_ratio,V = out$log_ratio_cov)
  #mash_data.L = mash_update_data(mash_data,ref=mash_ref)
  # canonical cov of mash
  U.c = cov_canonical(mash_data)
  # data-driven cov of mash
  m.1by1 = mash_1by1(mash_data)
  strong = get_significant_results(m.1by1)
  if(length(strong)>1){
    U.pca = cov_pca(mash_data,2,subset=strong)
    U.ed = cov_ed(mash_data, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data, c(U.c,U.ed), algorithm.version = 'Rcpp')
  }else{
    m = mash(mash_data, c(U.c), algorithm.version = 'Rcpp')
  }
  out_nonnull_agg[[simu]] = m

  save(out_nonnull_agg,file = '/scratch/midway2/dyxie/sc-cytokine/output/cd4_nonnull_pseudobulk.RData')





  ## directly use linear model

  vout = limma::voom(counts = thout$mat, design = X)
  fit_nonnull_voom = lash(vout,ref="1")
  out_nonnull_voom[[simu]] = fit_nonnull_voom

  save(out_nonnull_voom,file = '/scratch/midway2/dyxie/sc-cytokine/output/cd4_nonnull_voom.RData')




}





