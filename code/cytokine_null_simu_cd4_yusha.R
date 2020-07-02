library(Matrix)
################################ Data extraction #########################################
### load in the entire dataset (normalized and log-transformed) and sample annotation information
scdata <- readRDS("/scratch/midway2/dyxie/sc-cytokine/cyto_normalized.rds")
sample.info <- readRDS("/scratch/midway2/dyxie/sc-cytokine/cyto_samples.rds")

### reorder the columns by sequencing batch
samples <- as.character(unique(sample.info$sample))
samples <- c(samples[c(9,18,17,10,40)],samples[-c(9,18,17,10,40)])

### select data from a given cell type
scdata.CD4 <- t(scdata[sample.info$cell_type=="CD4_T_cells",])
sample.info.CD4 <- sample.info[sample.info$cell_type=="CD4_T_cells",]

### distribution of number of non-zero counts across genes
nonzero.counts <- apply(scdata.CD4,1,function(x){sum(x!=0)})
# hist(nonzero.counts, breaks="scott", xlab="Counts", main="Nonzero counts across genes in CD4")
# quantile(nonzero.counts, c(0.001,0.01,0.1))

# write a function to calculate Bhat and Shat
extract.data <- function(data,annot,trt){

  Bhat <- matrix(NA,nrow=nrow(data),ncol=length(trt))
  Shat <- matrix(NA,nrow=nrow(data),ncol=length(trt))
  rownames(Bhat) <- rownames(data)
  colnames(Bhat) <- trt
  rownames(Shat) <- rownames(data)
  colnames(Shat) <- trt

  overall.mean <- Matrix::rowMeans(data)

  for(i in 1:length(trt)){
    data.tmp <- data[, annot$sample==trt[i]]
    trt.mean <- Matrix::rowMeans(data.tmp)
    Bhat[,i] <- trt.mean
    Shat[,i] <- sqrt((apply(data.tmp, 1, var)+(trt.mean-overall.mean)^2)/ncol(data.tmp))
  }

  return(list(Bhat=Bhat, Shat=Shat))
}


rm(scdata)
######################### Apply mashr to the randomly permuted data #############################
library(mashr)

set.seed(1)

#pdf("mashr/simulations/Zscores_dist_combined.pdf", width=16, height=16)
#par(mfrow=c(4,4))

#simdata <- extract.data(scdata.CD4[nonzero.counts>=10, ], sample.info.CD4, samples)

# set up for mash
#data = mash_set_data(simdata$Bhat, simdata$Shat)
#data.L = mash_update_data(data, ref = 'mean')

#U.c = cov_canonical(data.L)
#m.1by1 = mash_1by1(data.L)
#strong = get_significant_results(m.1by1)
#U.pca = cov_pca(data.L,2,subset=strong)
#U.ed = cov_ed(data.L, U.pca, subset=strong)
#m = mash(data.L, c(U.c,U.ed), algorithm.version = 'Rcpp')


#hist(data.L$Bhat/data.L$Shat, breaks="scott", freq=FALSE, main="CD4 T cells: not permuted", xlab="Z scores", ylab="density", xlim=c(-10,5))

out_r250 = list()
out_r10 = list()
for (iter in 1:10){

  print(paste("running simulation: ", iter))

  ### remove the genes with very low counts across all treatments, use 10 for now
  # randomly permute the treatment labels to create null data
  idx = sample(1:ncol(scdata.CD4),ncol(scdata.CD4),replace=FALSE)
  scdata.perm <- scdata.CD4[nonzero.counts>=250, idx]

  # extract Bhat and Shat for the randomly permuted data
  simdata <- extract.data(scdata.perm, sample.info.CD4, samples)

  # set up for mash
  data = mash_set_data(simdata$Bhat, simdata$Shat)
  mash_data.L = mash_update_data(data, ref = 'mean')

  U.c = cov_canonical(mash_data.L)
  m.1by1 = mash_1by1(mash_data.L)
  strong = get_significant_results(m.1by1)
  if(length(strong)>1){
    U.pca = cov_pca(mash_data.L,2,subset=strong)
    U.ed = cov_ed(mash_data.L, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data.L, c(U.c,U.ed), algorithm.version = 'Rcpp',verbose = FALSE)
  }else{
    m = mash(mash_data.L, c(U.c), algorithm.version = 'Rcpp',verbose = FALSE)
  }
  out_r250[[iter]] = m
  print(paste("remove genes appearing in less than 250 cells: ",get_significant_results(m)))
  save(out_r250,file='output/cytokine_null_simu_cd4_r250_yusha.RData')


  ###################################
  ###################################
  scdata.perm <- scdata.CD4[nonzero.counts>=10, idx]

  # extract Bhat and Shat for the randomly permuted data
  simdata <- extract.data(scdata.perm, sample.info.CD4, samples)

  # set up for mash
  data = mash_set_data(simdata$Bhat, simdata$Shat)
  mash_data.L = mash_update_data(data, ref = 'mean')

  U.c = cov_canonical(mash_data.L)
  m.1by1 = mash_1by1(mash_data.L)
  strong = get_significant_results(m.1by1)
  if(length(strong)>1){
    U.pca = cov_pca(mash_data.L,2,subset=strong)
    U.ed = cov_ed(mash_data.L, U.pca, subset=strong)
    # fit mash model
    m = mash(mash_data.L, c(U.c,U.ed), algorithm.version = 'Rcpp',verbose = FALSE)
  }else{
    m = mash(mash_data.L, c(U.c), algorithm.version = 'Rcpp',verbose = FALSE)
  }
  out_r10[[iter]] = m
  print(paste("remove genes appearing in less than 10 cells: ",get_significant_results(m)))
  save(out_r10,file='output/cytokine_null_simu_cd4_r10_yusha.RData')

}
