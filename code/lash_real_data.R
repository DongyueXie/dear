library(limma)
library(mashr)
library(Matrix)
#'@param Y sample by feature matrix
#'@param groups treatment of each sample

lash = function(Y,groups){
  X = model.matrix(~as.factor(groups))
  X.lm = cbind(X[,1]-rowSums(X[,-1,drop=FALSE]),X[,-1])
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

## load dataset
samples = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_samples.rds')
datax = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_normalized.rds')
# cell_types = c('B_cells', 'CD4_T_cells', 'CD8_T_cells', 'NK_cells',
#                'Dendritic_cells','Ly6C+_Monocytes','Ly6C-_Monocytes','Neutrophils')

cell_types = c('CD4_T_cells', 'CD8_T_cells', 'NK_cells',
               'Dendritic_cells','Ly6C+_Monocytes','Ly6C-_Monocytes','Neutrophils')



# remove genes that have read counts in less than 'thresh' cells
thresh = 150

# for each cell type, run lash

out = list()

for(cell in cell_types){
  # extract data
  print(paste('running: ',cell))

  cell_idx = which(samples$cell_type == cell)
  cell_data = datax[cell_idx,]
  cell_samples = samples[cell_idx,]

  # remove genes
  rm_idx = which(colSums(cell_data!=0)<=thresh)
  if(length(rm_idx)>0){
    cell_data = cell_data[,-rm_idx]
  }
  out[[cell]] = lash(cell_data,cell_samples$sample)
  out[[cell]]$rm_gene_idx = rm_idx
  save(out,file = '/scratch/midway2/dyxie/sc-cytokine/output/lash_real_data_thresh150.RData')

}
