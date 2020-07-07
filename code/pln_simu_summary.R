library(mashr)
fdp = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    1-mean(dis.idx%in%true.idx)
  }
}


auc = function(pred,true.label){
  auc=pROC::roc(response = true.label, predictor = pred,direction = '<',levels = c(0,1))
  auc$auc
}

powr = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    sum(dis.idx%in%true.idx)/length(true.idx)
  }
}

summary_pln_simu = function(pln_simu,method.idx=c(1,2,3),null_names_simu = names(pln_simu)[1:2],
                            signal_names_simu=names(pln_simu)[3:4],skip_null=FALSE,non_null_mat_dim='full'){

  library(ggplot2)


  names_simu = names(pln_simu)

  if(!skip_null){

    cat("NULL case: compare with control and mean")
    cat('\n\n\n\n')

    null_summary = c()
    null_summary_sd = c()

    for(nn in null_names_simu){

      l = which(names_simu == nn)

      cat("Looking at",names_simu[l])
      methods = names(pln_simu[[l]])
      method = methods[method.idx]

      nsimu = length(pln_simu[[l]]$voom_limma)
      fdp_all = matrix(nrow=nsimu,ncol=length(method))
      n_fdis = matrix(nrow=nsimu,ncol=length(method))
      for(i in 1:nsimu){

        for(j in 1:length(method)){

          m = which(names(pln_simu[[l]])==method[j])

          fdp_all[i,j] = fdp(get_significant_results(pln_simu[[l]][[m]][[i]]$mash),0)
          n_fdis[i,j] = length(get_significant_results(pln_simu[[l]][[m]][[i]]$mash))

        }

      }
      colnames(fdp_all) = method
      colnames(n_fdis) = method
      rownames(n_fdis) = 1:nsimu
      null_summary = cbind(null_summary,c(apply(n_fdis,2,mean)))
      null_summary_sd = cbind(null_summary_sd,c(apply(n_fdis,2,sd)))
      print(knitr::kable(n_fdis,caption='number of false discoveries'))
      cat('\n\n\n\n')

    }
    colnames(null_summary) = c("control","mean")
    colnames(null_summary_sd) = c("control","mean")
    #print(knitr::kable(null_summary,caption='average #false discoveries'))
    #print(knitr::kable(null_summary_sd,caption='s.d. #false discoveries'))

    #cat('\n\n\n\n')


  }


  cat("Identifying significant genes")

  fdp_for_ggplot = data.frame()
  auc_for_ggplot = data.frame()
  power_for_ggplot = data.frame()

  for(nn in signal_names_simu){

    l = which(names_simu==nn)
    #cat("Looking at",names_simu[l])
    methods = names(pln_simu[[l]])
    method = methods[method.idx]

    nsimu = length(pln_simu[[l]]$voom_limma)

    fdp_signal = c()
    auc_signal = c()
    power_signal = c()

    for(i in 1:nsimu){

      which_null = 1*(rowSums(pln_simu[[l]]$non_null_matrix[[i]])==0)
      non_null_idx = which(rowSums(pln_simu[[l]]$non_null_matrix[[i]])!=0)

      for(j in 1:length(method)){

        m = which(names(pln_simu[[l]])==method[j])

        fdp_signal = rbind(fdp_signal,c(fdp(get_significant_results(pln_simu[[l]][[m]][[i]]$mash),non_null_idx) , method[j], unlist(strsplit(names_simu[l],split = '_'))[3]))
        auc_signal = rbind(auc_signal,c(auc(apply(pln_simu[[l]][[m]][[i]]$mash$result$lfsr,1,min),which_null) , method[j], unlist(strsplit(names_simu[l],split = '_'))[3]))
        power_signal = rbind(power_signal,c(powr(get_significant_results(pln_simu[[l]][[m]][[i]]$mash),non_null_idx), method[j], unlist(strsplit(names_simu[l],split = '_'))[3]))

      }

    }


    fdp_for_ggplot = rbind(fdp_for_ggplot,fdp_signal)
    auc_for_ggplot = rbind(auc_for_ggplot,auc_signal)
    power_for_ggplot = rbind(power_for_ggplot,power_signal)

    #colnames(fdp_signal) = method
    #colnames(auc_signal) = method
    #rownames(auc_signal) = 1:nsimu
    #colnames(power_signal) = method


    #boxplot(fdp_signal)
    #boxplot(auc_signal)
    #boxplot(power_signal)

    #print(knitr::kable(sort(round(apply(fdp_signal,2,mean),2)),col.names = 'lfsr level 0.05'))
    #print(knitr::kable(sort(round(apply(auc_signal,2,mean),2),decreasing = TRUE),col.names = 'auc'))
    #print(knitr::kable(sort(round(apply(power_signal,2,mean),2),decreasing = TRUE),col.names = 'power at lfsr level 0.05'))
    #cat('\n\n<!-- -->\n\n')

  }
  colnames(fdp_for_ggplot) = c("fdp","method","reference")
  colnames(auc_for_ggplot) = c("auc","method","reference")
  colnames(power_for_ggplot) = c("power","method","reference")


  fdp_for_ggplot$fdp = as.numeric(as.character(fdp_for_ggplot$fdp))
  auc_for_ggplot$auc = as.numeric(as.character(auc_for_ggplot$auc))
  power_for_ggplot$power = as.numeric(as.character(power_for_ggplot$power))


  p1 = ggplot(fdp_for_ggplot, aes(x=method, y=fdp, fill=reference)) +
    geom_boxplot() + theme(legend.position="none")

  p2 = ggplot(auc_for_ggplot, aes(x=method, y=auc, fill=reference)) +
    geom_boxplot() + theme(legend.position="none")

  p3 = ggplot(power_for_ggplot, aes(x=method, y=power, fill=reference)) +
    geom_boxplot()+ theme(legend.position="bottom")

  gridExtra::grid.arrange(p1,p2,p3,nrow=2)

  #browser()


  # identify  correct condition?

  cat("Identifying significant conditions")

  fdp_for_ggplot = data.frame()
  auc_for_ggplot = data.frame()
  power_for_ggplot = data.frame()

  for(nn in signal_names_simu){

    #cat("Looking at",names_simu[l])
    l = which(names_simu==nn)
    methods = names(pln_simu[[l]])
    method = methods[method.idx]

    nsimu = length(pln_simu[[l]]$voom_limma)

    fdp_for_ggplot = data.frame()
    auc_for_ggplot = data.frame()
    power_for_ggplot = data.frame()

    for(i in 1:nsimu){


      if(non_null_mat_dim=='full'){
        if(names_simu[l]=="signal_simu_control"){
        which_null = c(1-(pln_simu[[l]]$non_null_matrix[[i]])[,-1])
        non_null_idx = which((pln_simu[[l]]$non_null_matrix[[i]])[,-1] != 0)
      }else{
        which_null = c(1-pln_simu[[l]]$non_null_matrix[[i]])
        non_null_idx = which(pln_simu[[l]]$non_null_matrix[[i]]!=0)
      }
      }else{
        if(names_simu[l]=="signal_simu_control"){
          which_null = c(1-(pln_simu[[l]]$non_null_matrix[[i]]))
          non_null_idx = which((pln_simu[[l]]$non_null_matrix[[i]]) != 0)
        }else{
          which_null = c(1-cbind(rep(0,nrow(pln_simu[[l]]$non_null_matrix[[i]])),pln_simu[[l]]$non_null_matrix[[i]]))
          non_null_idx = which(cbind(rep(0,nrow(pln_simu[[l]]$non_null_matrix[[i]])),pln_simu[[l]]$non_null_matrix[[i]])!=0)
        }
      }



      #browser()

      for(j in 1:length(method)){

        m = which(names(pln_simu[[l]])==method[j])

        fdp_signal = rbind(fdp_signal,c(fdp(which(pln_simu[[l]][[m]][[i]]$mash$result$lfsr<0.05),non_null_idx), method[j], unlist(strsplit(names_simu[l],split = '_'))[3]) )
        auc_signal = rbind(auc_signal,c(auc(c(pln_simu[[l]][[m]][[i]]$mash$result$lfsr),which_null), method[j], unlist(strsplit(names_simu[l],split = '_'))[3]) )
        power_signal = rbind(power_signal,c(powr(which(pln_simu[[l]][[m]][[i]]$mash$result$lfsr<0.05),non_null_idx), method[j], unlist(strsplit(names_simu[l],split = '_'))[3]) )

      }

    }

    fdp_for_ggplot = rbind(fdp_for_ggplot,fdp_signal)
    auc_for_ggplot = rbind(auc_for_ggplot,auc_signal)
    power_for_ggplot = rbind(power_for_ggplot,power_signal)

    # colnames(fdp_signal) = method
    # colnames(auc_signal) = method
    # rownames(auc_signal) = 1:nsimu
    # colnames(power_signal) = method
    #
    #
    # print(knitr::kable(sort(round(apply(fdp_signal,2,mean),2)),col.names = 'lfsr level 0.05'))
    # print(knitr::kable(sort(round(apply(auc_signal,2,mean),2),decreasing = TRUE),col.names = 'auc'))
    # print(knitr::kable(sort(round(apply(power_signal,2,mean),2),decreasing = TRUE),col.names = 'power at lfsr level 0.05'))
    # cat('\n\n<!-- -->\n\n')

  }

  colnames(fdp_for_ggplot) = c("fdp","method","reference")
  colnames(auc_for_ggplot) = c("auc","method","reference")
  colnames(power_for_ggplot) = c("power","method","reference")


  fdp_for_ggplot$fdp = as.numeric(as.character(fdp_for_ggplot$fdp))
  auc_for_ggplot$auc = as.numeric(as.character(auc_for_ggplot$auc))
  power_for_ggplot$power = as.numeric(as.character(power_for_ggplot$power))

  p1 = ggplot(fdp_for_ggplot, aes(x=method, y=fdp, fill=reference)) +
    geom_boxplot() + theme(legend.position="none")

  p2 = ggplot(auc_for_ggplot, aes(x=method, y=auc, fill=reference)) +
    geom_boxplot() + theme(legend.position="none")

  p3 = ggplot(power_for_ggplot, aes(x=method, y=power, fill=reference)) +
    geom_boxplot()+ theme(legend.position="bottom")

  gridExtra::grid.arrange(p1,p2,p3,nrow=2)



}
