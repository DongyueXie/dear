source('code/scMash2.R')
load("data/bulk_cytokin_SI.RData")
mat = as.matrix(mat)



simu_study_bulkMash = function(mat,n_simu = 20,n_sub_gene = 2000,
                          n_condition = 10,n_replicates = 3,
                          signal_sd = 0.8,non_null_gene_prop=0.1,n_non_null_cond = 3,
                          printevery = TRUE, file.name = "pln_simu", data.driven.cov=TRUE,
                          seed=12345){

  set.seed(seed)
  n_sub_sample = n_condition*n_replicates
  condition = rep(1:n_condition,each=n_replicates)
  condition = as.factor(condition)
  designmat = model.matrix(~condition)[,-1]
  nC = ncol(designmat)



  null_simu_control = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),edgeR=list())
  null_simu_mean = list(voom_limma=list(),rle_voom_limma=list(),deseq=list())
  signal_simu_control = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),edgeR=list(),non_null_matrix=list(),coefmat=list())
  signal_simu_mean = list(voom_limma=list(),rle_voom_limma=list(),deseq=list(),non_null_matrix=list(),coefmat=list())

  for(simu in 1:n_simu){
    if(printevery){
      print(paste("running simulation:",simu))
    }

    mat1 = (mat)[,sample(1:ncol(mat),n_sub_sample)]
    mat1 = mat1[which(rowSums(mat1==0)==0),]
    mat1 = mat1[sample(1:nrow(mat1),n_sub_gene),]


    try({

      # null simulation

      null_simu_control$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization=NULL,ref=1,data.driven.cov=data.driven.cov)
      null_simu_control$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="RLE",ref=1,data.driven.cov=data.driven.cov)
      null_simu_control$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',normalization=NULL,ref=1,data.driven.cov=data.driven.cov)
      null_simu_control$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR',normalization=NULL)

      null_simu_mean$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization=NULL,ref='mean',data.driven.cov=data.driven.cov)
      null_simu_mean$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="RLE",ref='mean',data.driven.cov=data.driven.cov)
      null_simu_mean$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',normalization=NULL,ref='mean',data.driven.cov=data.driven.cov)

      ######################################################
      ######################################################

      # create signal

      coefmat = matrix(stats::rnorm(nC * n_sub_gene,0,signal_sd),
                       ncol = nC,
                       nrow = n_sub_gene)

      non_null_matrix = matrix(0,nrow = n_sub_gene,ncol = nC)
      non_null_matrix[,(nC-1):(nC-n_non_null_cond+1)] = 1

      nng = sample(1:n_sub_gene,round(n_sub_gene*non_null_gene_prop))
      non_null_matrix[-nng,] = 0

      coefmat = coefmat*non_null_matrix
      coefmat[,nC] = -rowSums(coefmat[,-nC])
      non_null_matrix = 1*(coefmat!=0)
      non_null_matrix = cbind(rep(0,n_sub_gene),non_null_matrix)

      thout = thin_diff(mat = as.matrix(mat1), design_fixed = designmat,  coef_fixed=coefmat)
      mat1 = thout$mat


      #signal_simu_mean$limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization=NULL)
      #signal_simu_mean$rle_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization="RLE")
      #signal_simu_mean$tmm_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref='mean',normalization="TMM")

      signal_simu_mean$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization=NULL,data.driven.cov=data.driven.cov)
      signal_simu_mean$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization="RLE",data.driven.cov=data.driven.cov)
      #signal_simu_mean$tmm_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref='mean',normalization="TMM")
      signal_simu_mean$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',ref='mean',data.driven.cov=data.driven.cov)
      signal_simu_mean$non_null_matrix[[simu]] = non_null_matrix
      signal_simu_mean$coefmat[[simu]] = coefmat


      #signal_simu_control$limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization=NULL)
      #signal_simu_control$rle_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization="RLE")
      #signal_simu_control$tmm_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',ref=1,normalization="TMM")

      signal_simu_control$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization=NULL,data.driven.cov=data.driven.cov)
      signal_simu_control$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization="RLE",data.driven.cov=data.driven.cov)
      #signal_simu_control$tmm_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',ref=1,normalization="TMM")
      signal_simu_control$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2',ref=1,data.driven.cov=data.driven.cov)
      #signal_simu_control$anova[[simu]] = simple_aov(mat1,condition,method = 'anova')
      #signal_simu_control$kw[[simu]] = simple_aov(mat1,condition,method = 'kw')
      signal_simu_control$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR')
      signal_simu_control$non_null_matrix[[simu]] = non_null_matrix
      signal_simu_control$coefmat[[simu]] = coefmat


    })

    simu_all = list(null_simu_control=null_simu_control,null_simu_mean=null_simu_mean,
                    signal_simu_control=signal_simu_control,signal_simu_mean=signal_simu_mean)
    save(simu_all,file = paste("/scratch/midway2/dyxie/sc-cytokine/output/",file.name,".RData",sep=''))

  }


}



# case 1: n_replicate = 3, use data-driven cov

simu_study_bulkMash(mat,file.name = "bulkMash_3rep")

# case 2: n_replicate = 3, do not use data-driven cov

simu_study_bulkMash(mat,file.name = "bulkMash_3rep_no_datadriven_cov",data.driven.cov = FALSE)

# case 3: n_replicate =10, use data-driven cov

simu_study_bulkMash(mat,file.name = "bulkMash_10rep",n_replicates = 10)

# case 4: n_replicate =10, do not use data-driven cov

simu_study_bulkMash(mat,file.name = "bulkMash_10repno_datadriven_cov",n_replicates = 10,data.driven.cov = FALSE)
