library(mashr)
library(tidyverse)


sub_cpm = readRDS('data/BR_log_cpm.RDS')


summary = sub_cpm %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column('name') %>%
        extract(name, c('organ', 'condition', 'mouse'), '(.*)_(.*)_(.*)') %>%
        group_by(condition) %>%
        summarize_at(vars(-organ, -mouse), list(mean = mean, se = ~(sd(.) / sqrt(n()))))

    rownames(summary) = summary$condition

    C_hat = select(summary, contains('_mean'))
    colnames(C_hat) = str_replace(colnames(C_hat), '_mean', '')
    C_hat = t(as.matrix(C_hat))

    S_hat = select(summary, contains('_se'))
    colnames(S_hat) = str_replace(colnames(S_hat), '_se', '')
    S_hat = t(as.matrix(S_hat))

    # remove genes with a zero in S_hat
    to_remove = unique(which(S_hat == 0, arr.ind = T)[,'row'])
    C_hat = C_hat[-to_remove,]
    S_hat = S_hat[-to_remove,]


    # set contrast
    data = mash_set_data(C_hat, S_hat)
    data.L = mash_update_data(data, 'Control')

    # setup cov
    U.c = cov_canonical(data.L)
    m.1by1 = mash_1by1(data.L)
    strong = get_significant_results(m.1by1,0.05)
    U.pca = cov_pca(data.L,5,subset=strong)
    U.ed = cov_ed(data.L, U.pca, subset=strong)

    lapply(U.c,function(X){eigen(X)$values})
    lapply(U.ed,function(X){eigen(X)$values})

    # mash!
    mashcontrast.model = mash(data.L, c(U.c,U.ed))
    save(mashcontrast.model)

    U.pca2 = cov_pca(data.L,2,subset=strong)
    U.ed2 = cov_ed(data.L, U.pca2, subset=strong)
    mashcontrast.model2 = mash(data.L, c(U.c,U.ed2))
    save(mashcontrast.model)

    mashcontrast.model3 = mash(data.L, c(U.c))
    save(mashcontrast.model)
