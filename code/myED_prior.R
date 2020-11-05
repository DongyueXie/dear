
library(mvtnorm)
library(mvnfast)

#'@title prepare data for fitting mash, after fitting myED
myMash_data = function(ed.res,data){

  N = nrow(data$Bhat)
  R = ncol(data$Bhat)
  K = length(ed.res$Ulist)

  sigma2 = ed.res$sigma2
  sl = length(sigma2)
  Ulist = ed.res$Ulist
  if(sl==1){
    Ulist = lapply(Ulist, function(z){z+diag(sigma2,R)})
  }else if(sl==R){
    Ulist = lapply(Ulist, function(z){z+diag(sigma2)})
  }else if(sl==K){
    for(k in 1:K){
      Ulist[[k]] = Ulist[[k]] + diag(sigma2[k],R)
    }
  }else if(sl==N){
    data$Shat = sqrt(data$Shat^2 + sigma2%*%t(rep(1,R)))
  }else{
    stop('invalid ed output')
  }

  return(list(data=data,Ulist=Ulist))

}

#'@title wrapper function for myED
#'@param data mash data object
myED_wrapper = function(data,V_init,sigma2_init,sigma2_type,subsets=NULL,w_init=NULL,fix_sigma2=FALSE,...){
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
      ed.res = myED(data$Bhat[subsets,],data$Shat[subsets,]^2,w_init,V_init,sigma2_init,sigma2_type,fix_sigma2,...)
    }else{
      if(!is.null(data$L)){
        ycovar = lapply(subsets, function(i) data$L %*% (data$Shat_orig[i,] * t(data$V * data$Shat_orig[i,])) %*% t(data$L) )
      }else{
        ycovar = lapply(subsets, function(i) data$Shat[i,] * t(data$V * data$Shat[i,]) )
      }
      ed.res = myED(data$Bhat[subsets,],ycovar,w_init,V_init,sigma2_init,sigma2_type,fix_sigma2,...)
    }
  }else{
    ed.res = myED(data$Bhat[subsets,],data$V[,,subsets],w_init,V_init,sigma2_init,sigma2_type,fix_sigma2,...)
  }

  ed.res
}

#'@param X data, a N by R matrix
#'@param S a matrix, each row is the variance of x_i; or an array of covariance matrices.
#'@param w_init,V_init,sigma2_init initial value of pi(a vector),V(a list),sigma^2(a vector or a scalar)
#'@param sigma2_type universal,universalD,mixture-specific,sample-specific.
#'@param fix_sigma2 If true, do not update sigma2.
#'@param alpha,beta shape and rate of inverse gamma prior of sigma2
#'@param xMean0 if False, also estimate mean of X
#'@param tol tolerance for convergence
#'@param maxiter maximum number of iterations to perform
myED = function(X,S,w_init,V_init,sigma2_init,sigma2_type='universal',fix_sigma2=FALSE,alpha=1,beta=NULL,
                xMean0=TRUE,tol=1e-5,maxiter=1e6,printevery = 5){

  N = nrow(X)
  R = ncol(X)
  K = length(V_init)
  w_curr = w_init
  V_curr = V_init
  sigma2_curr = sigma2_init

  if(is.null(beta)){
    if(is.matrix(S)){
      beta = R/sqrt(sum((X/S)^2)+N*R)
    }else{
      beta = R/sqrt(sum((X/t(apply(S,3,diag)))^2)+N*R)
    }
  }

  if(xMean0){
    m = rep(0,N)
  }

  #gl_temp = calcLogLik(X,S,V_curr,sigma2_curr,w_curr)
  #gamma_curr = gl_temp$gammas
  #log_liks = gl_temp$LogLik
  L_curr = mixtureLik(X,S,V_curr,sigma2_curr,m,sigma2_type)
  obj = calcobj(L_curr,w_curr,sigma2_curr,sigma2_type,alpha,beta)
  run_times = matrix(nrow=maxiter,ncol=4)
  for(iter in 1:maxiter){

    if(iter%%printevery==0){print(sprintf("running iteration %d, obj %f)",iter,obj[iter]))}

    t1 = Sys.time()
    gamma_curr = update_gamma(L_curr,w_curr)
    #print(colSums(gamma_curr))
    w_curr = update_w(gamma_curr)
    t2 = Sys.time()
    #print(w_curr)
    V_curr = update_V(X,S,V_curr,sigma2_curr,gamma_curr,sigma2_type)
    t3 = Sys.time()
    #print(V_curr)
    if(!fix_sigma2){
      sigma2_curr = update_randomeffect(X,S,V_curr,sigma2_curr,gamma_curr,alpha,beta,sigma2_type)
    }

    t4 = Sys.time()
    #gl_temp = calcLogLik(X,S,V_curr,sigma2_curr,w_curr)
    #gamma_curr = gl_temp$gammas
    #log_liks[iter+1] = gl_temp$LogLik

    L_curr = mixtureLik(X,S,V_curr,sigma2_curr,m,sigma2_type)
    obj[iter+1] = calcobj(L_curr,w_curr,sigma2_curr,sigma2_type,alpha,beta)
    t5 = Sys.time()

    rt = c(t1,t2,t3,t4,t5)
    run_times[iter,] = rt[-1] - rt[-5]
    #print(log_liks)
    if(abs(obj[iter+1]-obj[iter])<tol){
      print('done!')
      break
    }
  }
  if(iter==maxiter){
    warning('Algorithm did not converge.')
  }

  V_curr = lapply(V_curr,function(z){(z+t(z))/2})

  return(list(pi = w_curr,Ulist = V_curr,sigma2 = sigma2_curr,obj=obj,beta=beta,run_time = run_times))

}

is.monotone = function(x){all((x[-1]-x[-length(x)])>0)}

#'@title calculate objective function
calcobj = function(L,w,sigma2,sigma2_type,alpha,beta){

  if(alpha<=0|beta<=0){
    obj = sum(log(drop(L%*%w)))
  }else{
    N = nrow(L)

    if(sigma2_type=='universal'){
      obj = sum(log(drop(L%*%w)))+ N*dinvgamma(sigma2,alpha,beta,log.d = T)
    }
    if(sigma2_type=='universalD'){
      obj = sum(log(drop(L%*%w)))+ N*sum(dinvgamma(sigma2,alpha,beta,log.d = T))
    }
    if(sigma2_type=='mixture-specific'){
      d = dinvgamma(sigma2,alpha,beta,log.d = F)
      obj = sum(log(drop(t(t(L)*d)%*%w)))
    }
  }
  obj
}

#'@title density of inverse gamma distribution
dinvgamma = function(x,alpha,beta,log.d=T){
  if(log.d){
    d = alpha*log(beta)-log(gamma(alpha)) - (alpha+1)*log(x) - beta/x
  }else{
    d = beta^alpha/gamma(alpha)*(1/x)^(alpha+1)*exp(-beta/x)
  }
  d
}

#'@title update random effect function
#'
update_randomeffect = function(X,S,V,sigma2,gammas,alpha,beta,sigma2_type){

  if(sigma2_type=='universal'){
    sigma2_curr = update_sigma2_universal(X,S,V,sigma2,gammas,alpha,beta)
  }else if(sigma2_type=='universalD'){
    sigma2_curr = update_sigma2_universalD(X,S,V,sigma2,gammas,alpha,beta)
  }else if(sigma2_type=='mixture-specific'){
    sigma2_curr = update_sigma2_mixture(X,S,V,sigma2,gammas,alpha,beta)
  }else if(sigma2_type=='sample-specific'){
    sigma2_curr = update_sigma2_sample(X,S,V,sigma2,gammas)
  }else{
    stop('invalid initial value of sigma2')
  }
  sigma2_curr
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
      Si = S[,,i]
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
update_sigma2_mixture = function(X,S,V,sigma2,gammas,alpha,beta){
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
        Si = S[,,i]
      }
      inv_temp = sigma2[k]*chol2inv(chol(Vk+Si+diag(sigma2[k],R)))
      #hik = inv_temp%*%X[i,]
      #Hik = sigma2[k]*inv_temp%*%(Vk+Si)
      delta_k = delta_k + gammas[i,k]*(sum((inv_temp%*%X[i,])^2)+sum(inv_temp*(Vk+Si))+2*beta)
    }
    sigma2[k] = (delta_k)/((R+2*alpha+2)*gammas_k[k])
  }
  sigma2
}


#'@title update sigma2, universal, D
update_sigma2_universalD = function(X,S,V,sigma2,gammas,alpha,beta){

  #print(sigma2)

  N = nrow(X)
  R = ncol(X)
  K = length(V)

  Sigma2 = diag(sigma2)
  isSmat = is.matrix(S)
  Delta = 0
  for(i in 1:N){
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[,,i]
    }

    #hi = 0
    #Hi = 0
    for(k in 1:K){
      Vk = V[[k]]
      inv_temp = sigma2*chol2inv(chol(Vk+Si+Sigma2))
      #hik = (inv_temp)%*%X[i,]
      #Hik = (inv_temp)%*%(Vk+Si)
      Delta = Delta + gammas[i,k]*(((inv_temp)%*%X[i,])^2+diag(inv_temp%*%(Vk+Si)))
      #Delta = Delta + gammas[i,k]*(sum(hik^2)+sigma2*sum(inv_temp*(Vk+Si)))
      #hi = hi+gammas[i,k]*hik
      #Hi = Hi + gammas[i,k]*(tcrossprod(hik)+Hik)
    }
    #Hi = Hi - tcrossprod(hi)
    #Delta = Delta + crossprod(hi) + sum(diag(Hi))
  }
  c((Delta+2*N*beta)/(N*(1+2*alpha+2)))
}

#'@title update sigma2, universal, sigma^2*I
update_sigma2_universal = function(X,S,V,sigma2,gammas,alpha,beta){

  #print(sigma2)

  N = nrow(X)
  R = ncol(X)
  K = length(V)

  Sigma2 = diag(sigma2,R)
  isSmat = is.matrix(S)
  Delta = 0
  for(i in 1:N){
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[,,i]
    }

    #hi = 0
    #Hi = 0
    for(k in 1:K){
      Vk = V[[k]]
      inv_temp = sigma2*chol2inv(chol(Vk+Si+Sigma2))
      #hik = inv_temp%*%X[i,]
      #Hik = sigma2*inv_temp%*%(Vk+Si)
      Delta = Delta + gammas[i,k]*(sum((inv_temp%*%X[i,])^2)+sum(inv_temp*(Vk+Si)))
      #hi = hi+gammas[i,k]*hik
      #Hi = Hi + gammas[i,k]*(tcrossprod(hik)+Hik)
    }
    #Hi = Hi - tcrossprod(hi)
    #Delta = Delta + crossprod(hi) + sum(diag(Hi))
  }
  c((Delta+2*N*beta)/(N*(R+2*alpha+2)))
}

#'
#' #'@title cal obj function, use log sum exp trick
#' calcLogLik = function(X,S,V,sigma2,w){
#'
#'   N = nrow(X)
#'   R = ncol(X)
#'   K = length(V)
#'   gammas= matrix(nrow = N,ncol=K)
#'   isSmat = is.matrix(S)
#'   LogLik = 0
#'   for(i in 1:N){
#'
#'     if(isSmat){
#'       Si = diag(S[i,])
#'     }else{
#'       Si = S[,,i]
#'     }
#'
#'     if(length(sigma2)==K){
#'       Sigma_i = V
#'       for(k in 1:K){
#'         Sigma_i[[k]] = Sigma_i[[k]] + Si + diag(sigma2[k],R)
#'       }
#'     }else{
#'       if(length(sigma2)==1){
#'         sigma2_temp = diag(sigma2,R)
#'       }else if(length(sigma2)==R){
#'         sigma2_temp = diag(sigma2)
#'       }else if(length(sigma2)==N){
#'         sigma2_temp = diag(sigma2[i],R)
#'       }
#'       Sigma_i = lapply(V,function(z){z+Si+sigma2_temp})
#'     }
#'
#'     temp_i = calcLogLik_onesample(X[i,],Sigma_i,w)
#'     LogLik = LogLik + temp_i$LogLik
#'     gammas[i,] = temp_i$gamma
#'   }
#'   return(list(LogLik = LogLik,gammas=gammas))
#' }
#'
#' #'@return 1. log sum exps ; 2. gamma_i
#' calcLogLik_onesample = function(x,Sigma,w){
#'   R = length(x)
#'   l_temp = unlist(lapply(Sigma,function(z,x){-t(x)%*%solve(z)%*%x/2},x=x))
#'   w_temp = lapply(Sigma,function(z){log(w/sqrt((2*pi)^R*det(z)))})
#'   w_temp = diag(do.call(rbind,w_temp))
#'
#'   a = max((l_temp+w_temp))
#'   LogLik = a + log(sum(exp(w_temp+l_temp-a)))
#'   gamma_i = exp(l_temp+w_temp-a)/sum(exp(l_temp+w_temp-a))
#'   return(list(LogLik=LogLik,gamma=gamma_i))
#' }

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

update_V = function(X,S,V,sigma2,gammas,sigma2_type){
  N = nrow(X)
  R = ncol(X)
  K = length(V)

  gamma_k = colSums(gammas)
  #check if S is matrix
  isSmat = is.matrix(S)
  #check sigma2 mode

  if(sigma2_type=='universal'){
    sigma2_temp = diag(sigma2,R)
  }
  if(sigma2_type=='universalD'){
    sigma2_temp = diag(sigma2)
  }

  for(k in 1:K){
    #obtain b_ik and B_ik
    #print(gamma_k)

    #if(abs(gamma_k[k])>sqrt(.Machine$double.eps)){

      #b_ik = matrix()
      #B_ik = matrix()
      Vk = V[[k]]
      Deltak = 0
      for(i in 1:N){
        if(isSmat){
          Si = diag(S[i,])
        }else{
          Si = S[,,i]
        }

        if(sigma2_type=='mixture-specific'){
          sigma2_temp = diag(sigma2[k],R)
        }else if(sigma2_type=='sample-specific'){
          sigma2_temp = diag(sigma2[i],R)
        }


        inv_temp = Vk%*%chol2inv(chol(Vk+Si+sigma2_temp))
        #bik = inv_temp%*%X[i,]
        #Bik = inv_temp%*%(sigma2_temp+Si)
        #Deltak = Deltak + gammas[i,k]*(tcrossprod(bik)+Bik)
        Deltak = Deltak + gammas[i,k]*(tcrossprod(inv_temp%*%X[i,])+inv_temp%*%(sigma2_temp+Si))
        #Deltak = Deltak + gammas[i,k]*(inv_temp%*%(tcrossprod(X[i,])%*%inv_temp+sigma2_temp+Si))
      }
      V[[k]] = Deltak/gamma_k[k]

    #}

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
mixtureLik = function(X,S,V,sigma2,m,sigma2_type){
  N = nrow(X)
  R = ncol(X)
  K = length(V)
  L = matrix(nrow=N,ncol=K)

  #check if S is matrix
  isSmat = is.matrix(S)

  if(sigma2_type=='universal'){
    sigma2_temp = diag(sigma2,R)
  }
  if(sigma2_type=='universalD'){
    sigma2_temp = diag(sigma2)
  }

  for(i in 1:N){
    xi = X[i,]
    if(isSmat){
      Si = diag(S[i,])
    }else{
      Si = S[,,i]
    }
    for(k in 1:K){

      if(sigma2_type=='mixture-specific'){
        sigma2_temp = diag(sigma2[k],R)
      }else if(sigma2_type=='sample-specific'){
        sigma2_temp = diag(sigma2[i],R)
      }

      cov_temp = Si + V[[k]] + sigma2_temp
      L[i,k] = dmvn(xi,rep(m[i],R),sigma = cov_temp)
    }
  }
  L

}

# test
# X = simdata$Bhat
# S = simdata$Shat
# V = simdata$U.true
# mixtureLik(X,S,V,sigma2 = rep(1,100),m=rep(0,100))
