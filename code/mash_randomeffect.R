library(mashr)
library(mvtnorm)


n_conditions = function(data){ncol(data$Bhat)}

n_effects = function(data){nrow(data$Bhat)}

bovy_wrapper = function(data, Ulist_init, subset=NULL, ...){
  if(is.null(subset)){subset = 1:n_effects(data)}
  K = length(Ulist_init)
  R = n_conditions(data)
  pi_init = rep(1/K, K) # initial mix proportions
  D = ncol(data$V)
  if(all(data$V==diag(D))){
    ed.res = extreme_deconvolution(data$Bhat[subset,],
                                   data$Shat[subset,]^2,
                                   xamp = pi_init,
                                   xmean = matrix(0,nrow=K,ncol=R),
                                   xcovar = Ulist_init,
                                   fixmean = TRUE,
                                   ...)
  }else{
    if(!is.null(data$L)){
      ycovar = lapply(subset, function(i) data$L %*% (data$Shat_orig[i,] * t(data$V * data$Shat_orig[i,])) %*% t(data$L) )
    }else{
      ycovar = lapply(subset, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
    }
    ed.res = extreme_deconvolution(data$Bhat[subset,],
                                   ycovar,
                                   xamp = pi_init,
                                   xmean = matrix(0,nrow=K,ncol=R),
                                   xcovar = Ulist_init,
                                   fixmean = TRUE,
                                   ...)
  }
  return(list(pi = ed.res$xamp, Ulist = ed.res$xcovar, av_loglik = ed.res$avgloglikedata))
}

calc_loglikx = function(data,subset,pihat,Ulist,sigma){

  n = length(subset)
  loglik = 0
  for(i in subset){
    loglik = loglik + mixture_loglikx(data$Bhat[i,],data$Shat[i,],pihat,Ulist,sigma)
  }
  loglik
}

mixture_loglikx = function(x,shat,pihat,Ulist,sigma){

  K = length(pihat)
  p = length(x)
  lik = 0
  for(k in 1:K){
    #browser()
    lik = lik + pihat[k]*dmvnorm(x,sigma = Ulist[[k]]+diag(shat^2)+sigma^2*diag(p))
  }
  log(lik)
}



fdp = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    1-mean(dis.idx%in%true.idx)
  }
}


auc = function(pred,true.label){
  auc=pROC::roc(response = true.label, predictor = pred,direction = '<',levels = c(0,1))
  auc$auc
}

powr = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    sum(dis.idx%in%true.idx)/length(true.idx)
  }
}

mse = function(x,y){
  mean((x-y)^2)
}

summary_out = function(B,out = list(m.c=m.c,m.ed=m.ed,m.c.ed=m.c.ed,m.true=m.true),
                       alpha=0.05,criteria = 'lfsr'){

  # identify genes
  non_null_idx = which(rowSums(B)!=0)
  non_null_idx_c = which(B!=0)

  which_null = 1*(rowSums(B)==0)
  which_null_c = 1*(B==0)

  fdps = c()
  #aucs = c()
  powers = c()

  fdps_c = c()
  #aucs_c = c()
  powers_c = c()

  mses = c()
  log_liks = c()

  for(i in 1:length(out)){

    if(criteria=='lfsr'){
      fdps[i] = fdp(get_significant_results(out[[i]],thresh = alpha),non_null_idx)
      fdps_c[i] = fdp(which(out[[i]]$result$lfsr<alpha),non_null_idx_c)

      #aucs[i] = auc(c(apply(out[[i]]$result$lfsr,1,min)),which_null)
      #aucs_c[i] = auc(c(out[[i]]$result$lfsr),c(which_null_c))

      powers[i] = powr(get_significant_results(out[[i]],thresh = alpha),non_null_idx)
      powers_c[i] = powr(which(out[[i]]$result$lfsr<alpha),non_null_idx_c)

    }


    if(criteria=='lfdr'){
      fdps[i] = fdp(get_significant_results(out[[i]],thresh = alpha),non_null_idx)
      fdps_c[i] = fdp(which(out[[i]]$result$lfdr<alpha),non_null_idx_c)

      #aucs[i] = auc(c(apply(out[[i]]$result$lfdr,1,min)),which_null)
      #aucs_c[i] = auc(c(out[[i]]$result$lfdr),c(which_null_c))

      powers[i] = powr(get_significant_results(out[[i]],thresh = alpha),non_null_idx)
      powers_c[i] = powr(which(out[[i]]$result$lfdr<alpha),non_null_idx_c)

    }


    mses[i] = mse(B,out[[i]]$result$PosteriorMean)
    log_liks[i] = get_loglik(out[[i]])
  }

  find_genes = rbind(fdps,powers)
  rownames(find_genes) = c('fdp','power')
  colnames(find_genes) = names(out)

  find_cond = rbind(fdps_c,powers_c,mses,log_liks)
  rownames(find_cond) = c('fdp','power','mse','log_lik')
  colnames(find_cond) = names(out)

  return(list(find_genes=find_genes,find_cond=find_cond,mses=mses))

}

simu_study = function(simdata,subset){
  data = mash_set_data(simdata$Bhat,simdata$Shat)
  #m.1by1 = mash_1by1(data)
  #strong = get_significant_results(m.1by1)
  strong = subset
  U.c    = cov_canonical(data)
  U.pca  = cov_pca(data,5,strong)
  U.ed   = cov_ed(data,U.pca,strong)

  U.true = simdata$U.true
  m.c    = mash(data, U.c,verbose = F)
  m.ed   = mash(data, U.ed,verbose = F)
  m.c.ed = mash(data, c(U.c,U.ed),verbose = F)
  #m.c.ed.sparse = mash(data, c(U.c,U.ed.sparse),verbose = F)
  m.true = mash(data, U.true,verbose = F)
  out = list(m.true=m.true,m.c=m.c,m.ed=m.ed,m.c.ed=m.c.ed)
  out
}



simple_sims0 = function(nsamp = 100, err_sd = 0.01){
  ncond = 5
  b1 = rnorm(nsamp)
  B.1 = matrix(cbind(b1, b1, 0, 0, 0), nrow = nsamp, ncol = ncond)
  b2 = rnorm(nsamp)
  B.2 = matrix(cbind(0, 0, b2, b2, b2), nrow = nsamp, ncol = ncond)

  #B.id = matrix(rnorm(nsamp * ncond), nrow = nsamp, ncol = ncond)
  #B.zero = matrix(0, nrow = nsamp, ncol = ncond)

  B = rbind(B.1, B.2)
  Shat = matrix(err_sd, nrow = nrow(B), ncol = ncol(B))
  E = matrix(rnorm(length(Shat), mean = 0, sd = Shat), nrow = nrow(B), ncol = ncol(B))
  Bhat = B + E
  row_ids = paste0("effect_", 1:nrow(B))
  col_ids = paste0("condition_", 1:ncol(B))
  rownames(B) = row_ids
  colnames(B) = col_ids
  rownames(Bhat) = row_ids
  colnames(Bhat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  U = matrix(0,nrow=ncond,ncol=ncond)
  U2 = U
  U2[1:2,1:2] = 1
  U3 = U
  U3[3:5,3:5] = 1
  U.true = list(#U1 = matrix(0,nrow=ncond,ncol=ncond),
    U2=U2,
    U3=U3)
  #U4 = diag(ncond))
  return(list(B = B, Bhat = Bhat, Shat = Shat,U.true=U.true))
}


simple_sims01=function(nsamp = 500, S,ncond = 5){
  b1 = rnorm(nsamp)
  B.1 = matrix(cbind(b1, b1, 0, 0, 0), nrow = nsamp, ncol = ncond)
  b2 = rnorm(nsamp)
  B.2 = matrix(cbind(0, 0, b2, b2, b2), nrow = nsamp, ncol = ncond)

  #B.id = matrix(rnorm(nsamp * ncond), nrow = nsamp, ncol = ncond)
  #B.zero = matrix(0, nrow = nsamp, ncol = ncond)

  B = rbind(B.1, B.2)
  Shat = matrix(1, nrow = nrow(B), ncol = ncol(B))
  E = matrix(nrow = nrow(B), ncol = ncol(B))
  for(i in 1:nrow(B)){
    E[i,] = MASS::mvrnorm(1,mu=rep(0,ncol(B)),Sigma = S[,,i])
  }
  #E = matrix(rnorm(length(Shat), mean = 0, sd = Shat), nrow = nrow(B), ncol = ncol(B))
  Bhat = B + E
  row_ids = paste0("effect_", 1:nrow(B))
  col_ids = paste0("condition_", 1:ncol(B))
  rownames(B) = row_ids
  colnames(B) = col_ids
  rownames(Bhat) = row_ids
  colnames(Bhat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  U = matrix(0,nrow=ncond,ncol=ncond)
  U2 = U
  U2[1:2,1:2] = 1
  U3 = U
  U3[3:5,3:5] = 1
  U.true = list(#U1 = matrix(0,nrow=ncond,ncol=ncond),
    U2=U2,
    U3=U3)
  #U4 = diag(ncond))
  return(list(B = B, Bhat = Bhat, Shat = Shat,U.true=U.true))
}
