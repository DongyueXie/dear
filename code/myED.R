
library(mvtnorm)

#'@title prepare data for fitting mash, after fitting myED
myMash_data = function(ed.res,data){

  N = nrow(data$Bhat)
  R = ncol(data$Bhat)
  K = length(ed.res$Ulist)

  sigma2 = ed.res$sigma2
  sl = length(sigma2)
  if(sl==1){

    Ulist = ed.res$Ulist
    Ulist = lapply(Ulist, function(z){z+diag(sigma2,R)})

  }else if(sl==K){

    Ulist = ed.res$Ulist
    for(k in 1:K){
      Ulist[[k]] = Ulist[[k]] + diag(sigma2[k],R)
    }

  }else if(sl==N){

    data$Shat = sqrt(data$Shat^2 + sigma2%*%t(rep(1,R)))
    Ulist = ed.res$Ulist

  }else{
    stop('invalid ed output')
  }

  return(list(data=data,Ulist=Ulist))

}

#'@title wrapper function for myED
#'@param data mash data object
myED_wrapper = function(data,V_init,sigma2_init,subsets=NULL,w_init=NULL,xMean0=TRUE,...){
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

  if(all(data$V==diag(D))){
    ed.res = myED(data$Bhat[subsets,],data$Shat[subsets,]^2,w_init,V_init,sigma2_init,xMean0,...)
  }else{
    if(!is.null(data$L)){
      ycovar = lapply(subsets, function(i) data$L %*% (data$Shat_orig[i,] * t(data$V * data$Shat_orig[i,])) %*% t(data$L) )
    }else{
      ycovar = lapply(subsets, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
    }
    ed.res = myED(data$Bhat[subsets,],ycovar,w_init,V_init,sigma2_init,xMean0,...)
  }
  ed.res
}

#'@param X data, a N by R matrix
#'@param S a matrix, each row is the variance of x_i; or a list of covariance matrices.
#'@param w_init,V_init,sigma2_init initial value of pi(a vector),V(a list),sigma^2(a vector or a scalar)
#'@param xMean0 if False, also estimate mean of X
#'@param tol tolerance for convergence
#'@param maxiter maximum number of iterations to perform
myED = function(X,S,w_init,V_init,sigma2_init,xMean0=TRUE,tol=1e-5,maxiter=1e6,printevery = 5){

  N = nrow(X)
  K = length(V_init)
  w_curr = w_init
  V_curr = V_init
  sigma2_curr = sigma2_init

  if(xMean0){
    m = rep(0,N)
  }
  L_curr = mixtureLik(X,S,V_curr,sigma2_curr,m)
  log_liks = calcLogLik(L_curr,w_curr)
  for(iter in 1:maxiter){

    if(iter%%printevery==0){print(sprintf("running iteration %d, loglik %f)",iter,log_liks[iter]))}

    gamma_curr = update_gamma(L_curr,w_curr)
    w_curr = update_w(gamma_curr)
    V_curr = update_V(X,S,V_curr,sigma2_curr,gamma_curr)
    if(length(sigma2_init)==1){

      sigma2_curr = update_sigma2_universal(X,S,V_curr,sigma2_curr,gamma_curr)

    }else if(length(sigma2_init)==K){

      sigma2_curr = update_sigma2_mixture(X,S,V_curr,sigma2_curr,gamma_curr)

    }else if(length(sigma2_init)==N){

      sigma2_curr = update_sigma2_sample(X,S,V_curr,sigma2_curr,gamma_curr)

    }else{
      stop('invalid initial value of sigma2')
    }


    L_curr = mixtureLik(X,S,V_curr,sigma2_curr,m)
    log_liks[iter+1] = calcLogLik(L_curr,w_curr)
    if(abs(log_liks[iter+1]-log_liks[iter])<tol){
      print('done!')
      break
    }
  }
  if(iter==maxiter){
    warning('Algorithm did not converge.')
  }

  V_curr = lapply(V_curr,function(z){(z+t(z))/2})

  return(list(pi = w_curr,Ulist = V_curr,sigma2 = sigma2_curr,log_liks=log_liks))

}

#'@title update sigma2, sample specific
update_sigma2_sample = function(X,S,V,sigma2,gammas){
  N = nrow(X)
  R = ncol(X)
  K = length(V)

  isSmat = is.matrix(S)
  Delta = 0
  for(i in 1:N){
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[[i]]
    }
    xi = X[i,]

    hi = 0
    Hi = 0
    for(k in 1:K){
      Vk = V[[k]]
      inv_temp = solve(Vk+Si+diag(sigma2[i],R))
      hik = sigma2[i]*inv_temp%*%xi
      Hik = sigma2[i]*inv_temp%*%(Vk+Si)
      hi = hi+gammas[i,k]*hik
      Hi = Hi + gammas[i,k]*(tcrossprod(hik)+Hik)
    }
    Hi = Hi - tcrossprod(hi)
    sigma2[i] = (crossprod(hi)+sum(diag(Hi)))/R
  }
  sigma2
}

#'@title update sigma2, mixture specific
update_sigma2_mixture = function(X,S,V,sigma2,gammas){
  N = nrow(X)
  R = ncol(X)
  K = length(V)
  gammas_k = colSums(gammas)

  isSmat = is.matrix(S)
  for(k in 1:K){
    Vk = V[[k]]

    delta_k = 0
    for(i in 1:N){
      if(isSmat){
        Si = diag(S[i,])
      }else{
        Si = S[[i]]
      }
      inv_temp = solve(Vk+Si+diag(sigma2[k],R))
      hik = sigma2[k]*inv_temp%*%X[i,]
      Hik = sigma2[k]*inv_temp%*%(Vk+Si)
      delta_k = delta_k + gammas[i,k]*(crossprod(hik)+sum(diag(Hik)))
    }
    sigma2[k] = delta_k/(R*gammas_k[k])
  }
  sigma2
}


#'@title update sigma2, universal
update_sigma2_universal = function(X,S,V,sigma2,gammas){

  #print(sigma2)

  N = nrow(X)
  R = ncol(X)
  K = length(V)

  isSmat = is.matrix(S)
  Delta = 0
  for(i in 1:N){
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[[i]]
    }
    xi = X[i,]

    hi = 0
    Hi = 0
    for(k in 1:K){
      Vk = V[[k]]
      inv_temp = solve(Vk+Si+diag(sigma2,R))
      hik = sigma2*inv_temp%*%xi
      Hik = sigma2*inv_temp%*%(Vk+Si)
      hi = hi+gammas[i,k]*hik
      Hi = Hi + gammas[i,k]*(tcrossprod(hik)+Hik)
    }
    Hi = Hi - tcrossprod(hi)
    Delta = Delta + crossprod(hi) + sum(diag(Hi))
  }
  c(Delta/(N*R))
}

#'@title update gamma
#'@param L likelihood matrix
#'@param w weights
update_gamma = function(L,w){
  temp = t(t(L)*w)
  temp/rowSums(temp)
}
#'@title update weights
#'@param gammas posterior weights
update_w = function(gammas){
  colMeans(gammas)
}

#'@title update V

update_V = function(X,S,V,sigma2,gammas){
  N = nrow(X)
  R = ncol(X)
  K = length(V)

  gamma_k = colSums(gammas)
  #check if S is matrix
  isSmat = is.matrix(S)
  #check sigma2 mode
  if(length(sigma2)==N){
    sigma2 = sigma2%*%t(rep(1,K))
  }
  if(length(sigma2)==K){
    sigma2 = rep(1,N)%*%t(sigma2)
  }
  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=N,ncol=K)
  }

  for(k in 1:K){
    #obtain b_ik and B_ik
    b_ik = matrix()
    B_ik = matrix()
    Vk = V[[k]]
    Deltak = 0
    for(i in 1:N){
      if(isSmat){
        Si = diag(S[i,])
      }else{
        Si = S[[i]]
      }
      inv_temp = Vk%*%solve(Vk+Si+sigma2[i,k]*diag(R))
      bik = inv_temp%*%X[i,]
      Bik = inv_temp%*%(diag(sigma2[i,k],R)+Si)
      Deltak = Deltak + gammas[i,k]*(bik%*%t(bik)+Bik)
    }
    V[[k]] = Deltak/gamma_k[k]
  }
  V

}

#'@title calculate log-likelihood
#'@param L likelihood matrix, N*K
#'@param w mixture weights, length K vector
calcLogLik = function(L,w){
  #sum(log(rowSums(L%*%diag(w))))
  sum(log(drop(L%*%w)))
}

#'@title evaluate likelihood for each mixture component
#'@param X data, a N by R matrix
#'@param S a matrix, each row is the variance of x_i; or a list of covariance matrices.
#'@param V a list of length K, mixture component covariance matrix
#'@param sigma2 either a scalar or a length K vector or a length N vector
#'@param m the mean of each sample, length N vector
#'@return L, likelihood matrix
mixtureLik = function(X,S,V,sigma2,m){
  N = nrow(X)
  R = ncol(X)
  K = length(V)
  L = matrix(nrow=N,ncol=K)

  #check if S is matrix
  isSmat = is.matrix(S)
  #check sigma2 mode
  if(length(sigma2)==N){
    sigma2 = sigma2%*%t(rep(1,K))
  }
  if(length(sigma2)==K){
    sigma2 = rep(1,N)%*%t(sigma2)
  }
  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=N,ncol=K)
  }

  for(i in 1:N){
    xi = X[i,]
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[[i]]
    }
    for(k in 1:K){
      cov_temp = Si + V[[k]] + sigma2[i,k]*diag(R)
      L[i,k] = dmvnorm(xi,mean=rep(m[i],R),sigma = cov_temp)
    }
  }
  L

}

# test
# X = simdata$Bhat
# S = simdata$Shat
# V = simdata$U.true
# mixtureLik(X,S,V,sigma2 = rep(1,100),m=rep(0,100))
