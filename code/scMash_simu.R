source('code/scMash2.R')
# pbmc.data = Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
# pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 5, min.features = 200)
# pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
# pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# mat = pbmc@assays$RNA@counts
# load("data/pbmc3k_cellType.RData")
#
#
# mat = mat[,c(which(cell_type=='Naive CD4 T'),which(cell_type=='Memory CD4 T'),
#              which(cell_type=='CD14+ Mono'),which(cell_type=='B'))]



### cytokine data
samples = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_samples.rds')
datax = readRDS('/scratch/midway2/dyxie/sc-cytokine/data/cyto_counts.rds')
cell = "CD4_T_cells"
cell_idx = which(samples$cell_type == cell)
mat = t(datax[cell_idx,])
remove(datax)

mat = mat[which(rowSums(mat > 0) >= 100),]

# Null simulation
nsimu = 20
non_null_gene_prop = 0.1
non_null_trt_prop = 0.2
n_condition = 10
signal_sd = 0.8

n_sub_gene = 2000
n_sub_cell = 50*3*n_condition


set.seed(12345)

null_simu = list(limma=list(),
                 voom_limma=list(),
                 #deseq=list(),
                 limma_pb=list(),
                 voom_limma_pb=list(),
                 deseq_pb=list(),
                 anova=list())
signal_simu = list(limma=list(),
                   voom_limma=list(),
                   #deseq_a0=list(),
                   #limma_a1=list(),
                   #deseq_a1=list(),
                   limma_pb=list(),
                   voom_limma_pb=list(),
                   deseq_pb=list(),
                   #limma_pb_a1=list(),
                   #deseq_pb_a1=list(),
                   anova=list(),
                   non_null_matrix=list(),
                   coefmat=list())

for(simu in 1:nsimu){
  print(paste("running simulation:",simu))
  mat1 = (mat)[sample(1:nrow(mat),n_sub_gene),sample(1:ncol(mat),n_sub_cell)]
  condition = sample(1:n_condition,ncol(mat1),replace = TRUE)
  celltype = as.factor(condition)



try({
  null_simu$limma[[simu]] = mde_mash(mat1,celltype,test.method = 'limma',ref='mean',pseudoBulk = FALSE)
  null_simu$voom_limma[[simu]] = mde_mash(mat1,celltype,test.method = 'voom-limma',ref='mean',pseudoBulk = FALSE)

  null_simu$limma_pb[[simu]] = mde_mash(mat1,celltype,test.method = 'limma',ref='mean',pseudoBulk = T)
  null_simu$voom_limma_pb[[simu]] = mde_mash(mat1,celltype,test.method = 'voom-limma',ref='mean',pseudoBulk = T)
  null_simu$deseq_pb[[simu]] = mde_mash(mat1,celltype,test.method = 'DESeq2',ref='mean',pseudoBulk = T)

  null_simu$anova[[simu]] = simple_aov(mat1,celltype)


  # create signals
  designmat = model.matrix(~0+celltype)[,-1]
  G = nrow(mat1)
  nC = ncol(designmat)

  coefmat = matrix(stats::rnorm(nC * G,0,signal_sd),
                   ncol = nC,
                   nrow = G)

  non_null_matrix = matrix(0,nrow = G,ncol = nC)

  nng = sample(1:G,round(G*non_null_gene_prop))
  nnc = round(nC*non_null_trt_prop)
  for(j in nng){
    non_null_matrix[j,sample(1:nC,(nnc+sample(c(-1,0,1),1)))] = 1
  }

  coefmat = coefmat*non_null_matrix
  coefmat[,nC] = 0-rowSums(coefmat[,-nC])
  non_null_matrix = 1*(coefmat!=0)

  thout = thin_diff(mat = as.matrix(mat1), design_fixed = designmat,  coef_fixed  = coefmat)

  signal_simu$limma[[simu]] = mde_mash(thout$mat,celltype,test.method = 'limma',ref='mean',pseudoBulk = FALSE,alpha=0)
  signal_simu$voom_limma[[simu]] = mde_mash(thout$mat,celltype,test.method = 'voom-limma',ref='mean',pseudoBulk = FALSE,alpha=0)
  #signal_simu$limma_a1[[simu]] = mde_mash(thout$mat,celltype,test.method = 'limma',ref='mean',pseudoBulk = FALSE,alpha=1)
  #signal_simu$deseq_a0[[simu]] = mde_mash(thout$mat,celltype,test.method = 'DESeq2',ref='mean',pseudoBulk = FALSE,alpha=0)
  #signal_simu$deseq_a1[[simu]] = mde_mash(thout$mat,celltype,test.method = 'DESeq2',ref='mean',pseudoBulk = FALSE,alpha=1)

  signal_simu$limma_pb[[simu]] = mde_mash(thout$mat,celltype,test.method = 'limma',ref='mean',pseudoBulk = T,alpha=0)
  signal_simu$voom_limma_pb[[simu]] = mde_mash(thout$mat,celltype,test.method = 'voom-limma',ref='mean',pseudoBulk = T,alpha=0)
  #signal_simu$limma_pb_a1[[simu]] = mde_mash(thout$mat,celltype,test.method = 'limma',ref='mean',pseudoBulk = T,alpha=1)
  signal_simu$deseq_pb[[simu]] = mde_mash(thout$mat,celltype,test.method = 'DESeq2',ref='mean',pseudoBulk = T,alpha=0)
  #signal_simu$deseq_pb_a1[[simu]] = mde_mash(thout$mat,celltype,test.method = 'DESeq2',ref='mean',pseudoBulk = T,alpha=1)

  signal_simu$anova[[simu]] = simple_aov(thout$mat,celltype)
  signal_simu$non_null_matrix[[simu]] = non_null_matrix
  signal_simu$coefmat[[simu]] = coefmat
})


  save(null_simu,file='/scratch/midway2/dyxie/sc-cytokine/output/scMash_cytokine_cd4_null_50CellPerRep.RData')
  save(signal_simu,file='/scratch/midway2/dyxie/sc-cytokine/output/scMash_cytokine_cd4_signal_50CellPerRep.RData')
}
