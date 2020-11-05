
library(daarem)
#'@title wrapper function for myED
#'@param data mash data object
myED_daarem_wrapper = function(data,V_init,sigma2_init,subsets=NULL,w_init=NULL,fix_sigma2=FALSE,...){
  N = nrow(data$Bhat)
  R = ncol(data$Bhat)
  if(is.null(subsets)){subsets = 1:N}
  K = length(V_init)
  if(is.null(w_init)){
    w_init  = rep(1/K,K) # initial mix proportions
  }
  if(length(sigma2_init)==N){
    sigma2_init = sigma2_init[subsets]
  }
  D = ncol(data$V)

  if(prod(dim(data$V))==(D^2)){
    if(all(data$V==diag(D))){
      ed.res = myED_daarem(data$Bhat[subsets,],data$Shat[subsets,]^2,w_init,V_init,sigma2_init,fix_sigma2,...)
    }else{
      if(!is.null(data$L)){
        ycovar = lapply(subsets, function(i) data$L %*% (data$Shat_orig[i,] * t(data$V * data$Shat_orig[i,])) %*% t(data$L) )
      }else{
        ycovar = lapply(subsets, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
      }
      ed.res = myED_daarem(data$Bhat[subsets,],ycovar,w_init,V_init,sigma2_init,fix_sigma2,...)
    }
  }else{
    ed.res = myED_daarem(data$Bhat[subsets,],data$V[,,subsets],w_init,V_init,sigma2_init,fix_sigma2,...)
  }

  ed.res
}


# params: w, sigma2, V's
myED_daarem = function(X,S,w_init,V_init,sigma2_init,fix_sigma2=FALSE,alpha=1,beta=NULL,
                       xMean0=TRUE,tol=1e-5,maxiter=1e6,order=10){

  N = nrow(X)
  R = ncol(X)
  K = length(V_init)
  params_init = c(w_init,sigma2_init,list2vector(V_init))

  if(is.null(beta)){
    if(is.matrix(S)){
      beta = R/sqrt(sum((X/S)^2)+N*R)
    }else{
      beta = R/sqrt(sum((X/t(apply(S,3,diag)))^2)+N*R)
    }
  }

  t1 = Sys.time()
  out = suppressWarnings(
    daarem(params_init,myED_updates_daarem,myED_obj_daarem,X,S,N,K,R,length(sigma2_init),alpha,beta,fix_sigma2,
           control = list(maxiter = maxiter,order = order,tol = tol,
                          mon.tol = 0.05,kappa = 20,alpha = 1.2)))

  params_out = out$par

  t2 = Sys.time()
  return(list(pi = params_out[1:K],
              Ulist = vector2list(params_out[(K+length(sigma2_init)+1):length(params_out)],R),
              sigma2 = params_out[(K+1):(K+length(sigma2_init))],
              log_liks=out$objfn.track,
              run_time = t2-t1))

}

myED_obj_daarem = function(params,X,S,N,K,R,len_sigma2,alpha,beta,fix_sigma2){
  L = mixtureLik(X,S,vector2list(params[(K+len_sigma2+1):length(params)],R),params[(K+1):(K+len_sigma2)],rep(0,N))
  calcLogLik(L,params[1:K])
}

myED_updates_daarem = function(params,X,S,N,K,R,len_sigma2,alpha,beta,fix_sigma2){

  # input x_t, out put x_t+1

  w = params[1:K]
  sigma2 = params[(K+1):(K+len_sigma2)]
  V = vector2list(params[(K+len_sigma2+1):length(params)],R)
  L = mixtureLik(X,S,V,sigma2,rep(0,N))
  # EM
  gamma_curr = update_gamma(L,w)
  w_curr = update_w(gamma_curr)
  V_curr = update_V(X,S,V,sigma2,gamma_curr)
  if(!fix_sigma2){
    if(len_sigma2==1){

      sigma2_curr = update_sigma2_universal(X,S,V_curr,sigma2,gamma_curr,alpha,beta)

    }else if(len_sigma2==R){

      sigma2_curr = update_sigma2_universalD(X,S,V_curr,sigma2,gamma_curr,alpha,beta)

    }else if(len_sigma2==K){

      sigma2_curr = update_sigma2_mixture(X,S,V_curr,sigma2,gamma_curr,alpha,beta)

    }else if(len_sigma2==N){

      sigma2_curr = update_sigma2_sample(X,S,V_curr,sigma2,gamma_curr)

    }else{
      stop('invalid initial value of sigma2')
    }
  }else{
    sigma2_curr = sigma2
  }


  return(c(w_curr,sigma2_curr,list2vector(V_curr)))

}

vector2list = function(x,R){
  nmat = length(x)/R^2
  X = list()
  for(i in 1:nmat){
    X[[i]] = matrix(x[((i-1)*R^2+1):(i*R^2)],ncol=R)
  }
  X
}

list2vector = function(X){
  unlist(lapply(X,c))
}

