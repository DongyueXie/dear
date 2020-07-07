library(seqgendiff)
library(limma)
library(DESeq2)
library(edgeR)
library(Matrix)
library(mashr)


#'@param mat gene by cell count matrix
#'@param condition indicate each cell's condition, for example, treatment or cell type.
#'@param normalization RLE or TMM or NULL for running limma
#'@param test.method DEseq2, limma, anova for performing DE on pseudo-bulk data
#'@param ref reference level of condition, default to be 'mean'; must be in one of condition
#'@param num.reps number of replicates to create for each condition


mde_mash = function(mat,condition,
                    normalization = "RLE", test.method = 'DESeq2',
                    ref='mean',pseudoBulk=FALSE,num.reps=3,verbose=FALSE,data.driven.cov=T){

  G = nrow(mat)
  N = ncol(mat)
  out = list()
  out$input = list(mat=mat,condition=condition,normalization=normalization,
                   test.method=test.method,
                   ref=ref,pseudoBulk=pseudoBulk,num.reps=num.reps,data.driven.cov=data.driven.cov)
  V = NULL

  condition = as.factor(condition)
  if(ref!='mean'){
    condition = relevel(condition,ref)
  }
  cond_level = levels(condition)
  n_cond = length(cond_level)
  cov_of_interest = 1:n_cond



  # create pseudo-bulk data

  if(pseudoBulk){
    Y = matrix(nrow=G,ncol=num.reps*n_cond)
    cond = c()
    for(i in 1:n_cond){
      cond_idx = which(condition==cond_level[i])
      reps_idx = sample(1:num.reps,length(cond_idx),replace = TRUE)
      Y[,((i-1)*num.reps+1):(i*num.reps)] = pseudo_bulk(mat[,cond_idx],reps_idx)
      cond = c(cond,rep(cond_level[i],num.reps))
    }
  }else{
    Y = mat
    cond = condition
  }

  X = model.matrix(~0+cond)

  if(test.method=='voom_limma'){

    if(!is.null(normalization)){
      Y = DGEList(Y)
      Y = calcNormFactors(Y,method = normalization)
    }

    dds = voom(Y,X)
    lmout = lmFit(object = dds, design = X)
    eout  = eBayes(lmout)
    Bhat = lmout$coefficients[,cov_of_interest,drop=FALSE]
    Shat = lmout$stdev.unscaled[,cov_of_interest,drop=FALSE] * sqrt(eout$s2.post)
    colnames(Bhat) = cond_level
    colnames(Shat) = cond_level

  }

  if(test.method=='voom_limma_V'){

    if(!is.null(normalization)){
      Y = DGEList(Y)
      Y = calcNormFactors(Y,method = normalization)
    }
    X = model.matrix(~cond)
    dds = voom(Y,design = X)
    lmout = lmFit(object = dds, design = X)
    eout  = eBayes(lmout)
    Bhat = lmout$coefficients[,-1,drop=FALSE]
    Shat = 1

    V = array(dim=c(n_cond-1,n_cond-1,G))
    for(i in 1:G){
      V[,,i] = (eout$s2.post[i] * solve(t(X)%*%diag(dds$weights[1,])%*%X))[-1,-1]
    }

    #colnames(Bhat) = cond_level

  }


  if(test.method=='limma'){


    if(!is.null(normalization)){
      if(normalization=='RLE'){
        Y = t(t(Y)/calcNormFactors(Y,method = 'RLE'))
      }
      if(normalization=='TMM'){
        Y = t(t(Y)/calcNormFactors(Y,method = 'TMM'))
      }
    }

    lmout = lmFit(object = log(1e6*(Y+0.5)%*%diag(c(1/colSums(Y)))), design = X)
    eout  = eBayes(lmout)
    Bhat = lmout$coefficients[,cov_of_interest,drop=FALSE]
    Shat = lmout$stdev.unscaled[,cov_of_interest,drop=FALSE] * sqrt(eout$s2.post)
    colnames(Bhat) = cond_level
    colnames(Shat) = cond_level


  }



  if(test.method=='DESeq2'){
    cond = cbind(cond)
    dds = DESeqDataSetFromMatrix(countData = Y,colData = cond,design = X)
    dds = DESeq(dds,quiet=!verbose)
    res_names = resultsNames(dds)[cov_of_interest]
    Bhat = matrix(nrow = G,ncol=n_cond)
    Shat = matrix(nrow = G,ncol=n_cond)
    colnames(Bhat) = cond_level
    colnames(Shat) = cond_level
    for(i in 1:n_cond){
      res = results(dds,name=res_names[i])
      Bhat[,i] = res$log2FoldChange
      Shat[,i] = res$lfcSE
    }
  }



  Bhat[is.na(Bhat)] = 0
  Shat[is.na(Shat)] = Inf

  out$mash = mash_wrapper(Bhat,Shat,ref=ref,V=V,verbose=verbose,data.driven.cov=data.driven.cov)

  return(out)

}


mash_wrapper = function(Bhat,Shat,ref='mean',V=NULL,verbose=FALSE,data.driven.cov=FALSE){

  if(is.null(V)){

    data = mash_set_data(Bhat, Shat)
    data.L = mash_update_data(data, ref = ref)
    U.c = cov_canonical(data.L)
    if(data.driven.cov){
      m.1by1 = mash_1by1(data.L)
      strong = get_significant_results(m.1by1)
      #browser()
      if(length(strong)>2){
        U.pca = cov_pca(data.L,2,subset=strong)
        U.ed = cov_ed(data.L, U.pca, subset=strong)
      }else{
        U.ed = NULL
      }
    }else{
      U.ed = NULL
    }

    m = mash(data.L, c(U.c,U.ed), algorithm.version = 'Rcpp',verbose = verbose)
    m$input = list(Bhat=Bhat,Shat=Shat,ref=ref)
    m

  }else{

    data = mash_set_data(Bhat, V=V)

    U.c = cov_canonical(data)
    m = mash(data, c(U.c), algorithm.version = 'R',verbose = verbose)

  }


}



pseudo_bulk = function(mat,reps){
  reps = as.factor(reps)
  n_reps = length(unique(reps))
  n_gene = nrow(mat)
  reps_name = levels(reps)
  Xa = matrix(nrow=n_gene,ncol=n_reps)
  for(i in 1:n_reps){
    group_idx = which(reps==reps_name[i])
    Xa[,i] = rowSums(mat[,group_idx,drop=FALSE])
  }
  colnames(Xa) = reps_name
  Xa
}



BH = function(p,alpha=0.05){
  n=length(p)
  i=rank(p)
  idx = which(p<=(i/n*alpha))
  if(length(idx)==0){
    NULL
  }else{
    i0= max(i[idx])
    rej.idx = which(i<=i0)
    rej.idx
  }
}

simple_aov = function(mat,condition,method='anova',normalization='RLE'){



  if(method=='edgeR'){

    y = DGEList(counts = mat)
    y <- calcNormFactors(y)
    X = model.matrix(~condition)
    y <- estimateDisp(y,design = X)
    fit <- glmQLFit(y, X)
    results <- glmQLFTest(fit)
    p_val = results$table$PValue
    input = list(method=method,normalization='TMM')
  }else{


    if(!is.null(normalization)){
      if(normalization=='RLE'){
        mat = t(t(mat)/calcNormFactors(mat,method = 'RLE'))
      }
      if(normalization=='TMM'){
        mat = t(t(mat)/calcNormFactors(mat,method = 'TMM'))
      }
    }


    p_val = c()

    if(method=='anova'){

      mat = log(mat+0.5)

      for (i in 1:nrow(mat)) {
        fit = aov(y~.,data.frame(y=as.numeric(mat[i,]),x=condition))
        p_val[i] = summary(fit)[[1]][["Pr(>F)"]][1]
      }

    }

    if(method=='kw'){

      for (i in 1:nrow(mat)) {
        fit = kruskal.test(y~.,data.frame(y=as.numeric(mat[i,]),x=condition))
        p_val[i] = fit$p.value
      }

    }

    input = list(method=method,normalization=normalization)

  }





  return(list(input=input,p_val=p_val,sig_idx = BH(p_val)))
}


