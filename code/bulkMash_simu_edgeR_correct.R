source('code/scMash2.R')

# datax = read.csv('/project2/mstephens/chevrier-stephens/data/bulk_cytokine_filtered.csv',row.names = 1)
# samples = colnames(datax)
# samples = matrix(unlist(strsplit(samples,split = '_')),ncol=3,byrow = TRUE)
#
#
# sum(datax==0)/prod(dim(datax))
#
# ## use data from organ SI
#
# mat = datax[,which(samples[,1]=='SI')]
# mat = mat[-which(rowSums(mat>0)==0),]
load("data/bulk_cytokin_SI.RData")
mat = as.matrix(mat)

#mat = mat[-which(rowSums(mat==0)>0),]
# treatment = samples[which(samples[,1]=='SI'),2]

# rm(datax)

# method: anova, voom-limma, rle-voom-limma, tmm-voom-limma, deseq2

load('/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_null_dense2.RData')
#load('/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_signal_mean_ref_dense2.RData')
load('/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_signal_control_ref_dense2.RData')
## permute data, check null setting
### in every repetition, randomly draw 2000 genes and 30 samples, randomly assign one of ten labels to each sample.

n_simu = 20
n_condition = 10
n_sub_gene = 2000
n_sub_sample = n_condition*3

non_null_gene_prop = 0.1
non_null_trt_prop = 0.2
signal_sd = 0.8

condition = rep(1:n_condition,each=3)
condition = as.factor(condition)

set.seed(12345)

for(simu in 1:n_simu){
  print(paste("running simulation:",simu))



  try({

    # null simulation

    # null_simu$limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',normalization=NULL)
    # null_simu$rle_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',normalization="RLE")
    # null_simu$tmm_limma[[simu]] = mde_mash(mat1,condition,test.method = 'limma',normalization="TMM")
    #
    # null_simu$voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization=NULL)
    # null_simu$rle_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="RLE")
    # null_simu$tmm_voom_limma[[simu]] = mde_mash(mat1,condition,test.method = 'voom_limma',normalization="TMM")
    # null_simu$deseq[[simu]] = mde_mash(mat1,condition,test.method = 'DESeq2')
    #
    # null_simu$anova[[simu]] = simple_aov(mat1,condition,method = 'anova')
    # null_simu$kw[[simu]] = simple_aov(mat1,condition,method = 'kw')
    mat1 = null_simu$limma[[simu]]$input$mat
    null_simu$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR')


    ######################################################
    ######################################################

    # create signal

    mat1 = signal_simu_control$limma[[simu]]$input$mat


    # signal_simu_mean$limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref='mean',normalization=NULL)
    # signal_simu_mean$rle_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref='mean',normalization="RLE")
    # signal_simu_mean$tmm_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref='mean',normalization="TMM")
    #
    # signal_simu_mean$voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref='mean',normalization=NULL)
    # signal_simu_mean$rle_voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref='mean',normalization="RLE")
    # signal_simu_mean$tmm_voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref='mean',normalization="TMM")
    # signal_simu_mean$deseq[[simu]] = mde_mash(thout$mat,condition,test.method = 'DESeq2',ref='mean')
    # signal_simu_mean$non_null_matrix[[simu]] = non_null_matrix
    # signal_simu_mean$coefmat[[simu]] = coefmat
    #
    #
    # signal_simu_control$limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref=1,normalization=NULL)
    # signal_simu_control$rle_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref=1,normalization="RLE")
    # signal_simu_control$tmm_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'limma',ref=1,normalization="TMM")
    #
    # signal_simu_control$voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref=1,normalization=NULL)
    # signal_simu_control$rle_voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref=1,normalization="RLE")
    # signal_simu_control$tmm_voom_limma[[simu]] = mde_mash(thout$mat,condition,test.method = 'voom_limma',ref=1,normalization="TMM")
    # signal_simu_control$deseq[[simu]] = mde_mash(thout$mat,condition,test.method = 'DESeq2',ref=1)
    # signal_simu_control$anova[[simu]] = simple_aov(thout$mat,condition,method = 'anova')
    # signal_simu_control$kw[[simu]] = simple_aov(thout$mat,condition,method = 'kw')
    signal_simu_control$edgeR[[simu]] = simple_aov(mat1,condition,method = 'edgeR')
    # signal_simu_control$non_null_matrix[[simu]] = non_null_matrix
    # signal_simu_control$coefmat[[simu]] = coefmat


  })



  # no dense: only remove genes appear in less than n_condition
  # dense: removes genes do not appear in more than 3 sample
  # dense2: dense and add limma, edgeR, voom(Y) dropped X

  save(null_simu,file='/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_null_dense2.RData')
  #save(signal_simu_mean,file='/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_signal_mean_ref_dense2.RData')
  save(signal_simu_control,file='/scratch/midway2/dyxie/sc-cytokine/output/bulkMash_cytokine_SI_signal_control_ref_dense2.RData')


}



## add signals, check power, fdr, auc, etc.
