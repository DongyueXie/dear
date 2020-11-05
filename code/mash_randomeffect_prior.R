source('code/myED_prior.R')
source('code/mash_randomeffect.R')

## diagonal s, 0.1^2

set.seed(12345)
simdata = simple_sims0(500,err_sd = 0.1)

data = mash_set_data(simdata$Bhat,simdata$Shat)
#m.1by1 = mash_1by1(data)
#strong = get_significant_results(m.1by1)
strong = 1:nrow(simdata$B)
U.c    = cov_canonical(data)
U.pca  = cov_pca(data,5,strong)
#ed.out = bovy_wrapper(data,U.pca,strong)

myed.universal = myED_wrapper(data,U.pca,0.1,sigma2_type = 'universal',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universal,file='output/myed_universal01.RData')
myed.universal.prior = myED_wrapper(data,U.pca,0.1,maxiter=2000,sigma2_type = 'universal',printevery = 100,beta=NULL,alpha=1)
save(myed.universal.prior,file='output/myed_universal_prior01.RData')

myed.universalD = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universalD,file='output/myed_universalD01.RData')
myed.universalD.prior = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.universalD.prior,file='output/myed_universalD_prior01.RData')

myed.mixture = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.mixture,file='output/myed_mixture01.RData')
myed.mixture.prior = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.mixture.prior,file='output/myed_mixture_prior01.RData')


## diagonal s, 0.5^2

set.seed(12345)
simdata = simple_sims0(500,err_sd = 0.5)

data = mash_set_data(simdata$Bhat,simdata$Shat)
#m.1by1 = mash_1by1(data)
#strong = get_significant_results(m.1by1)
strong = 1:nrow(simdata$B)
U.c    = cov_canonical(data)
U.pca  = cov_pca(data,5,strong)
#ed.out = bovy_wrapper(data,U.pca,strong)

myed.universal = myED_wrapper(data,U.pca,0.1,sigma2_type = 'universal',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universal,file='output/myed_universal05.RData')
myed.universal.prior = myED_wrapper(data,U.pca,0.1,maxiter=2000,sigma2_type = 'universal',printevery = 100,beta=NULL,alpha=1)
save(myed.universal.prior,file='output/myed_universal_prior05.RData')

myed.universalD = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universalD,file='output/myed_universalD05.RData')
myed.universalD.prior = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.universalD.prior,file='output/myed_universalD_prior05.RData')

myed.mixture = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.mixture,file='output/myed_mixture05.RData')
myed.mixture.prior = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.mixture.prior,file='output/myed_mixture_prior05.RData')


## diagonal s, 0.5^2,change sample size

set.seed(12345)
simdata = simple_sims0(1000,err_sd = 0.5)

data = mash_set_data(simdata$Bhat,simdata$Shat)
#m.1by1 = mash_1by1(data)
#strong = get_significant_results(m.1by1)
strong = 1:nrow(simdata$B)
U.c    = cov_canonical(data)
U.pca  = cov_pca(data,5,strong)
#ed.out = bovy_wrapper(data,U.pca,strong)

myed.universal = myED_wrapper(data,U.pca,0.1,sigma2_type = 'universal',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universal,file='output/myed_universal05n1000.RData')
myed.universalD = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universalD,file='output/myed_universalD05n1000.RData')
myed.mixture = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.mixture,file='output/myed_mixture05n1000.RData')



myed.universal.prior = myED_wrapper(data,U.pca,0.1,maxiter=2000,sigma2_type = 'universal',printevery = 100,beta=NULL,alpha=1)
save(myed.universal.prior,file='output/myed_universal_prior05n1000.RData')

myed.universalD.prior = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.universalD.prior,file='output/myed_universalD_prior05n1000.RData')

myed.mixture.prior = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.mixture.prior,file='output/myed_mixture_prior05n1000.RData')


## nondiagonal s,



set.seed(12345)
n = 1000
S = rWishart(n,10,diag(5)/20)
simdata = simple_sims01(500,S)
data = mash_set_data(simdata$Bhat,simdata$Shat,V=S)
strong = 1:nrow(simdata$B)
U.c    = cov_canonical(data)
U.pca  = cov_pca(data,5,strong)


myed.universal = myED_wrapper(data,U.pca,0.1,sigma2_type = 'universal',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universal,file='output/myed_universalS.RData')
myed.universalD = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.universalD,file='output/myed_universalDS.RData')
myed.mixture = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=0,alpha=-1)
save(myed.mixture,file='output/myed_mixtureS.RData')



myed.universal.prior = myED_wrapper(data,U.pca,0.1,maxiter=2000,sigma2_type = 'universal',printevery = 100,beta=NULL,alpha=1)
save(myed.universal.prior,file='output/myed_universal_priorS.RData')

myed.universalD.prior = myED_wrapper(data,U.pca,rep(0.1,5),sigma2_type = 'universalD',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.universalD.prior,file='output/myed_universalD_priorS.RData')

myed.mixture.prior = myED_wrapper(data,U.pca,rep(0.1,length(U.pca)),sigma2_type = 'mixture-specific',maxiter=2000,printevery = 100,beta=NULL,alpha=1)
save(myed.mixture.prior,file='output/myed_mixture_priorS.RData')
