source('code/scMash2.R')

#'@param X design matrix, N by n_condition
#'@param B coef matrix, G by n_condition
#'@param Sigma variance, either a vector of length G or a G*G matrix
#'@param s library size of a sample is rpois(s)*G
generate_pln = function(X,B,Sigma,s=50){

  N = nrow(X)
  G = nrow(B)
  if(is.null(dim(Sigma))){
    E = matrix(rnorm(N*G,0,sqrt(Sigma)),nrow=N,byrow=TRUE)
  }else{
    E = mvnfast::rmvn(N,rep(0,G),Sigma)
  }
  S = rpois(N,s)*G
  M = exp(X%*%t(B))
  M = M/rowSums(M)
  M = exp(log(M*S)+E)
  Y = matrix(rpois(G*N,M),nrow=N)
  t(Y)
}


simu_study_pln = function(n_simu = 20,n_sub_gene = 2000,
                          n_condition = 10,n_replicates = 3,
                          Sigma = rep(0.1,n_sub_gene),s=50,
                          signal_sd = 0.8,non_null_gene_prop=0.1,n_non_null_cond = 3,
                          printevery = TRUE, file.name = "pln_simu",
                          seed=12345){

  set.seed(seed)
  condition = rep(1:n_condition,each=n_replicates)
  condition = as.factor(condition)
  X = model.matrix(~condition)
  n_sub_sample = n_condition*n_replicates


  null_simu_control = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),edgeR=list())
  null_simu_mean = list(voom_limma=list(),rle_voom_limma=list(),deseq=list())
  signal_simu_control = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),edgeR=list(),non_null_matrix=list(),coefmat=list())
  signal_simu_mean = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),non_null_matrix=list(),coefmat=list())


  B_null = matrix(0,nrow=n_sub_gene,ncol=n_condition)
  B_null[,1] = 1

  for(simu in 1:n_simu){
    if(printevery){
      print(paste("running simulation:",simu))
    }


    try({

      # null simulation


      mat1 = generate_pln(X,B_null,Sigma,s)

      null_simu_control$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization=NULL,ref=1)
      null_simu_control$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="RLE",ref=1)
      null_simu_control$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',normalization=NULL,ref=1)
      null_simu_control$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR',normalization=NULL)

      null_simu_mean$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization=NULL,ref='mean')
      null_simu_mean$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="RLE",ref='mean')
      null_simu_mean$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',normalization=NULL,ref='mean')

      ######################################################
      ######################################################

      # create signal

      B = matrix(stats::rnorm(n_sub_gene * n_condition,0,signal_sd),
                 ncol = n_condition,
                 nrow = n_sub_gene)

      non_null_matrix = matrix(0,nrow = n_sub_gene,ncol = n_condition)
      non_null_matrix[,(n_condition-1):(n_condition-n_non_null_cond+1)] = 1

      nng = sample(1:n_sub_gene,round(n_sub_gene*non_null_gene_prop))
      non_null_matrix[-nng,] = 0

      B = B*non_null_matrix
      B[,n_condition] = -rowSums(B)
      non_null_matrix = (1*(B!=0))
      B[,1] = 1


      mat1 = generate_pln(X,B,Sigma,s)



      #signal_simu_mean$limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization=NULL)
      #signal_simu_mean$rle_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization="RLE")
      #signal_simu_mean$tmm_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization="TMM")

      signal_simu_mean$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization=NULL)
      signal_simu_mean$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization="RLE")
      #signal_simu_mean$tmm_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization="TMM")
      signal_simu_mean$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',ref='mean')
      signal_simu_mean$non_null_matrix[[simu]] = non_null_matrix
      signal_simu_mean$coefmat[[simu]] = B


      #signal_simu_control$limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization=NULL)
      #signal_simu_control$rle_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization="RLE")
      #signal_simu_control$tmm_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization="TMM")

      signal_simu_control$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization=NULL)
      signal_simu_control$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization="RLE")
      #signal_simu_control$tmm_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization="TMM")
      signal_simu_control$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',ref=1)
      #signal_simu_control$anova[[simu]] = simple_aov(mat1,condition,method = 'anova')
      #signal_simu_control$kw[[simu]] = simple_aov(mat1,condition,method = 'kw')
      signal_simu_control$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR')
      signal_simu_control$non_null_matrix[[simu]] = non_null_matrix
      signal_simu_control$coefmat[[simu]] = B


    })

    pln_simu = list(null_simu_control=null_simu_control,null_simu_mean=null_simu_mean,
                    signal_simu_control=signal_simu_control,signal_simu_mean=signal_simu_mean)
    save(pln_simu,file = paste("output/",file.name,".RData",sep=''))

  }


}


# case1: diagnal Sigma, all equal

simu_study_pln(file.name = "pln_simu_diagSigma01")


# case2: diagnal Sigma, from data,
Sigma = diag(sigma(myPLN))
Sigma[which(Sigma>0.6)] = 0.3
simu_study_pln(file.name = "pln_simu_diagSigma_plnFit",Sigma=Sigma)

# case3: dense Sigma, from data,
Sigma = (sigma(myPLN))
simu_study_pln(file.name = "pln_simu_denseSigma_plnFit",Sigma=Sigma)


