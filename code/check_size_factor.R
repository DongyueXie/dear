pseudo_bulk = function(X,groups){
  n_group = length(unique(groups))
  n_gene = ncol(X)
  group_name = unique(groups)
  Xa = matrix(nrow=n_group,ncol=n_gene)
  for(i in 1:n_group){
    group_idx = which(groups==group_name[i])
    Xa[i,] = colSums(X[group_idx,])
  }
  rownames(Xa) = group_name
  Xa
}



samples = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_samples.rds')
datax = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_counts.rds')
# cell_types = c('B_cells', 'CD4_T_cells', 'CD8_T_cells', 'NK_cells',
#                'Dendritic_cells','Ly6C+_Monocytes','Ly6C-_Monocytes','Neutrophils')

cell_types = c('B_cells','CD4_T_cells', 'CD8_T_cells', 'NK_cells',
               'Dendritic_cells','Ly6C+_Monocytes','Ly6C-_Monocytes','Neutrophils')

size_factors = c()

for(cell in cell_types){
  cell_idx = which(samples$cell_type == cell)
  cell_data = datax[cell_idx,]
  cell_data = cell_data[,-which(colSums(cell_data)==0)]
  pseudo_bulk_data = pseudo_bulk(cell_data,samples$sample[cell_idx])


  mean_ge = apply(pseudo_bulk_data,2,mean)

  ###
  gre = log(diag(1/rowSums(pseudo_bulk_data))%*%(pseudo_bulk_data+0.5) )
  gre = gre[,order(mean_ge,decreasing = TRUE)]
  X = model.matrix(~rowSums(pseudo_bulk_data))
  lmout = limma::lmFit(t(gre),design = X)
  bin_idx = ceiling(seq_along(1:nrow(lmout))/(nrow(lmout)/20))
  plot_data = data.frame(y=lmout$coefficients[,2],grp=bin_idx)
  boxplot(y~grp,plot_data,xlab='bins',ylab='coefficients',main=paste(cell))



  ###
  size_factors = cbind(size_factors,rowSums(pseudo_bulk_data))
}
colnames(size_factors) = cell_types
boxplot(size_factors,ylab='Size Factor')


