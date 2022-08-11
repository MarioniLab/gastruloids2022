suppressPackageStartupMessages(library(reticulate))
use_condaenv("/nfs/research/marioni/Leah/software/conda_envs/deconvolution")
source_python('/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/perturbation_simulations/DA_functions.py')

io <- list()

##########
## npcs ##
##########
io$npcs <- list()

## d3 ##
#io$npcs[["d3"]] <- list()
#io$npcs[["d3"]][["All"]] = list("mesodermal"=3, "neural"=2)
#io$npcs[["d3"]][["exp4_d3"]] = list("mesodermal"=3, "neural"=2)

## d3.5 ##
#io$npcs[["d3.5"]] <- list()
#io$npcs[["d3.5"]][["All"]] = list("mesodermal"=3, "neural"=2)
#io$npcs[["d3.5"]][["exp5_d3.5_d4"]] = list("mesodermal"=3, "neural"=2)

## d4 ##
io$npcs[["d4"]] <- list()
io$npcs[["d4"]][["All"]] = list("mesodermal"=3, "neural"=3)
io$npcs[["d4"]][["exp3B_d4_d4.5"]] = list("mesodermal"=3, "neural"=3)
io$npcs[["d4"]][["exp5_d3.5_d4"]] = list("mesodermal"=3, "neural"=3)

## d4.5 ##
io$npcs[["d4.5"]] <- list()
io$npcs[["d4.5"]][["All"]] = list("mesodermal"=3, "neural"=2)
io$npcs[["d4.5"]][["exp3B_d4_d4.5"]] = list("mesodermal"=3, "neural"=2)
io$npcs[["d4.5"]][["exp6_d4.5_d5"]] = list("mesodermal"=3, "neural"=2)

## d5 ##
#io$npcs[["d5"]] <- list()
#io$npcs[["d5"]][["All"]] = list("mesodermal"=3, "neural"=2, "intermediate"=0)
#io$npcs[["d5"]][["exp1_d5"]] = list("mesodermal"=3, "neural"=2, "intermediate"=0)
#io$npcs[["d5"]][["exp6_d4.5_d5"]] = list("mesodermal"=3, "neural"=2, "intermediate"=0)

#############
## tps_ref ##
#############
io$tps_ref <- list()
io$tps_ref[["d3"]] <- c("d3", "d3.5")
io$tps_ref[["d3.5"]] <- c("d3", "d3.5", "d4")
io$tps_ref[["d4"]] <- c("d3.5", "d4", "d4.5")
io$tps_ref[["d4.5"]] <- c("d4", "d4.5", "d5")
io$tps_ref[["d5"]] <- c("d4.5", "d5")

io$gastr_type_subset <- list()
io$gastr_type_subset[["d3"]] <- c("mesodermal", "neural")
io$gastr_type_subset[["d3.5"]] <- c("mesodermal", "neural")
io$gastr_type_subset[["d4"]] <- c("mesodermal", "neural")
io$gastr_type_subset[["d4.5"]] <- c("mesodermal", "neural")
io$gastr_type_subset[["d5"]] <- c("mesodermal", "neural", "intermediate")

sample_gampois <- function(p_ct, norgs, cells_per_samp=10000, IoD, IoD_contr=0) {
    if (is.null(IoD)) {
        stop('IoD must be provided when sim_type==gampois')
    }
    if (IoD_contr == 0) {
        org_contr <- rep((cells_per_samp/norgs), norgs)
    } else if (IoD_contr == 1) {
        org_contr <- rpois(norgs, lambda=cells_per_samp/norgs)
    } else {
        # Note: IoD doesn't translate 1:1 because it's not scale invariant
        org_contr <- rgamma(norgs, shape=(cells_per_samp/norgs)/IoD_contr, scale=IoD_contr)
    }
    return(colSums(lapply(org_contr, function(n) unlist(lapply(p_ct, function(p) rpois(1, lambda=rgamma(1, shape=p*n/IoD, scale=IoD))))) %>% do.call(rbind, .)))
}

sample_poisson <- function(p_ct, norgs, IoD_contr=0, cells_per_samp=10000) {
    if (IoD_contr == 0) {
        org_contr <- rep((cells_per_samp/norgs), norgs)
    } else if (IoD_contr == 1) {
        org_contr <- rpois(norgs, lambda=cells_per_samp/norgs)
    } else {
        # Note: IoD doesn't translate 1:1 because it's not scale invariant
        org_contr <- rgamma(norgs, shape=(cells_per_samp/norgs)/IoD_contr, scale=IoD_contr)
    }
    return(colSums(lapply(org_contr, function(n) unlist(lapply(p_ct, function(p) rpois(1, lambda=p*n)))) %>% do.call(rbind, .)))
}

sample_multinom <- function(p_ct, norgs, IoD_contr=0, cells_per_samp=10000) {
    if (IoD_contr == 0) {
        org_contr <- rep((cells_per_samp/norgs), norgs)
    } else if (IoD_contr == 1) {
        org_contr <- rpois(norgs, lambda=cells_per_samp/norgs)
    } else {
        # Note: IoD doesn't translate 1:1 because it's not scale invariant
        org_contr <- rgamma(norgs, shape=(cells_per_samp/norgs)/IoD_contr, scale=IoD_contr)
    }
    return(rowSums(lapply(org_contr, function(n) rmultinom(1, size=n, prob=p_ct)) %>% do.call(cbind, .)))
}

sample_full <- function(to.sample, norgs, fc_ct, IoD_contr=0, cells_per_samp=10000) {
    if (IoD_contr == 0) {
        org_contr <- rep((cells_per_samp/norgs), norgs)
    } else if (IoD_contr == 1) {
        org_contr <- rpois(norgs, lambda=cells_per_samp/norgs)
    } else {
        # Note: IoD doesn't translate 1:1 because it's not scale invariant
        org_contr <- rgamma(norgs, shape=(cells_per_samp/norgs)/IoD_contr, scale=IoD_contr)
    }
    p_idx <- sample(1:nrow(to.sample), size=norgs, replace=TRUE)
    return(lapply(1:norgs, function(g) {
        p_ct <- to.sample[p_idx[g],]
        return(round(org_contr[g]*(p_ct*fc_ct/sum(p_ct*fc_ct))))
    } ) %>% do.call(rbind, .) %>% colSums)
}

do_experiment_general <- function(ns, fc_ct, p_ct, norgs, IoD, IoD_contr, sim_type, modality, cells_per_samp=10000) {
    if ((sim_type=="poisson") & (modality=="unimodal")) {
        trial_wt <- lapply(1:ns, function(i) sample_poisson(p_ct=p_ct, cells_per_samp=cells_per_samp, norgs=norgs, IoD_contr=IoD_contr)) %>% as.data.table 
        trial_perturb <- lapply(1:ns, function(i) sample_poisson(p_ct=(p_ct*fc_ct/sum(p_ct*fc_ct)), cells_per_samp=cells_per_samp, norgs=norgs, IoD_contr=IoD_contr)) %>% as.data.table 
        trial <- cbind(trial_wt, trial_perturb)
    }
    if ((sim_type=="gampois") & (modality=="unimodal")) {
        trial_wt <- lapply(1:ns, function(i) sample_gampois(p_ct=p_ct, cells_per_samp=cells_per_samp, norgs=norgs, IoD=IoD, IoD_contr=IoD_contr)) %>% as.data.table 
        trial_wt <- lapply(1:ns, function(i) sample_gampois(p_ct=(p_ct*fc_ct/sum(p_ct*fc_ct)), cells_per_samp=cells_per_samp, norgs=norgs, IoD=IoD, IoD_contr=IoD_contr)) %>% as.data.table 
        trial <- cbind(trial_wt, trial_perturb)
    }
    if ((modality!="unimodal")) stop("Only unimodal implemented right now")
    trial$celltype <- paste0("Celltype", 1:length(p_ct))
    trial <- melt(trial, variable.name="sample", value.name="ncells", id.vars="celltype")
    trial[1:(nrow(trial)/2), condition := "wt"]
    trial[((nrow(trial)/2)+1):nrow(trial), condition := "perturb"]
    if (sim_type=="gampois") trial$IoD <- IoD else trial$IoD <- NA
    if (!(is.null(IoD_contr))) trial$IoD_contr <- IoD_contr else trial$IoD_contr <- NA
    trial$sim_type <- sim_type
    trial$modality <- modality
    trial$norgs <- norgs
    trial$ns <- ns
    trial$nct <- length(p_ct)
    trial$cells_per_samp <- cells_per_samp
    return(trial)
}

do_experiment_gastruloid <- function(to.sample, fc_ct, ns, ng, IoD_contr, IoD=NULL, sim_type, cells_per_samp=10000) {
    if ((sim_type=="poisson")) {
        p_ct <- colMeans(to.sample)
        trial_wt <- lapply(1:ns, function(i) sample_poisson(p_ct=p_ct, cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial_perturb <- lapply(1:ns, function(i) sample_poisson(p_ct=(p_ct*fc_ct/sum(p_ct*fc_ct)), cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial <- cbind(trial_wt, trial_perturb)
    }
    if ((sim_type=="gampois")) {
        p_ct <- colMeans(to.sample)
        trial_wt <- lapply(1:ns, function(i) sample_gampois(p_ct=p_ct, cells_per_samp=cells_per_samp, norgs=ng, IoD, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial_perturb <- lapply(1:ns, function(i) sample_gampois(p_ct=(p_ct*fc_ct/sum(p_ct*fc_ct)), cells_per_samp=cells_per_samp, norgs=ng, IoD, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial <- cbind(trial_wt, trial_perturb)
    }
    if ((sim_type=="multinom")) {
        p_ct <- colMeans(to.sample)
        trial_wt <- lapply(1:ns, function(i) sample_multinom(p_ct=p_ct, cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial_perturb <- lapply(1:ns, function(i) sample_multinom(p_ct=(p_ct*fc_ct/sum(p_ct*fc_ct)), cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial <- cbind(trial_wt, trial_perturb)
    }
    if ((sim_type=="full")) {
        #p_ct <- colMeans(to.sample)
        trial_wt <- lapply(1:ns, function(i) sample_full(to.sample=to.sample, fc_ct=rep(1, ncol(to.sample)), cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial_perturb <- lapply(1:ns, function(i) sample_full(to.sample=to.sample, fc_ct=fc_ct, cells_per_samp=cells_per_samp, norgs=ng, IoD_contr=IoD_contr)) %>% do.call(cbind, .)
        trial <- cbind(trial_wt, trial_perturb)
    }
    cts <- rownames(trial)
    trial <- as.data.table(trial)
    trial$celltype <- cts
    trial <- melt(trial, variable.name="sample", value.name="ncells", id.vars="celltype")
    trial[1:(nrow(trial)/2), condition := "wt"]
    trial[((nrow(trial)/2)+1):nrow(trial), condition := "perturb"]
    if (sim_type=="gampois") trial$IoD <- IoD else trial$IoD <- NA
    if (!(is.null(IoD_contr))) trial$IoD_contr <- IoD_contr else trial$IoD_contr <- NA
    trial$sim_type <- sim_type
    trial$ng <- ng
    trial$ns <- ns
    trial$cells_per_samp <- cells_per_samp
    return(trial)
}

do_poisson_DA <- function(trial) {
    df <- group_by(trial, sample) %>%
        mutate(N_s=sum(ncells)) %>%
        ungroup() %>%
        group_by(celltype) %>%
        do(model=glm("ncells ~ condition + offset(log(N_s))", data=., family=stats::poisson()))
    res_df <- t(sapply(df$model, function(x) summary(x)$coefficients[nrow(summary(x)$coefficients),]))
    colnames(res_df) <- c("logFC","Std. Error", "z value",    "Pval" )
    louvain.res <- cbind(df, res_df) %>%
        mutate(FDR=p.adjust(Pval, method = "BH"))
    louvain.res$model <- NULL
    return(louvain.res)
}

do_NB_DA <- function(trial, norm.method='TMM') {
    design_df <- as_tibble(trial[,c("sample", "condition")]) %>%
        distinct()

    to.model <- dcast(trial, celltype ~ sample, value.var="ncells")
    rns <- to.model$celltype
    to.model$celltype <- NULL
    to.model <- as.matrix(to.model)
    rownames(to.model) <- rns

    if(norm.method %in% c("TMM")){
        dge <- DGEList(counts=to.model,
                       lib.size=colSums(to.model))
        dge <- calcNormFactors(dge, method="TMM")
    } else if(norm.method %in% c("logMS")){
        dge <- DGEList(counts=to.model,
                       lib.size=colSums(to.model))
    } else {
        stop('Unknown norm.method')
    }

    model <- model.matrix(formula('~ condition'), data=design_df)
    rownames(model) <- design_df$sample

    dge <- estimateDisp(dge, model)
    fit <- glmQLFit(dge, model, robust=TRUE)
    n.coef <- ncol(model)
    louvain.res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
    louvain.res$celltype <- rownames(louvain.res)
    louvain.res <- as.data.table(louvain.res)
    louvain.res <- setnames(louvain.res, "PValue", "Pval")
    return(louvain.res)
}

do_Deconv_DA <- function(trial, io, pca) {
    
    y_err <- list()
    ncts <- dim(pca$X)[2]
    for (ct_idx in 1:ncts) {

        y_err[[ct_idx]] <- list()
        samps <- as.character(unique(trial$sample))
        for (samp in samps) {

            toc <- matrix(trial[sample == samp]$ncells)
            rownames(toc) <- trial[sample == samp]$celltype
            toc <- as.matrix(toc[pca$ct_ordering,])

            if (ct_idx == 1) {
                model_out <- run_model(as.matrix(toc[2:ncts,]),
                                       pca$X[,2:ncts,],
                                       pca$offset[2:ncts,]
                                      )
            } else if (ct_idx < ncts) {
                model_out <- run_model(as.matrix(toc[c(1:(ct_idx-1),(ct_idx+1):ncts),]),
                                       pca$X[,c(1:(ct_idx-1),(ct_idx+1):ncts),],
                                       pca$offset[c(1:(ct_idx-1),(ct_idx+1):ncts),]
                                      )
            } else {
                model_out <- run_model(as.matrix(toc[1:(ct_idx-1),]),
                                       pca$X[,1:(ct_idx-1),],
                                       pca$offset[1:(ct_idx-1),]
                                      )
            }
            y_pred <- rowSums((exp(apply((model_out[[2]][, rep(1, dim(pca$X)[2]), ] * pca$X), MARGIN=c(2,3), sum) + pca$offset)) * abind(lapply(1:ncts, function(j) t(model_out[[1]])), along=1))
            y_err[[ct_idx]][[samp]] <- toc - y_pred
        }

        y_err[[ct_idx]] <- Reduce(cbind, y_err[[ct_idx]])
        colnames(y_err[[ct_idx]]) <- samps

    }

    design <- unique(trial[,c("sample", "condition")]) %>% as.data.frame
    rownames(design) <- design$sample
    design$sample <- NULL
    t.test.list <- lapply(pca$ct_ordering, function(ct) t.test(y_err[[which(pca$ct_ordering==ct)]][ct,] ~ design[colnames(y_err[[which(pca$ct_ordering==ct)]]),"condition"], alternative = "two.sided", paired=FALSE, var.equal=TRUE))
    names(t.test.list) <- pca$ct_ordering
    pvals <- lapply(t.test.list, function(t) t$p.value) %>% unlist
    mean_perturb <- lapply(t.test.list, function(t) t$estimate[["mean in group perturb"]]) %>% unlist
    mean_wt <- lapply(t.test.list, function(t) t$estimate[["mean in group wt"]]) %>% unlist
    out <- data.table(celltype = pca$ct_ordering, Pval = pvals, FDR = p.adjust(pvals, method="BH"), mean_perturb = mean_perturb, mean_wt = mean_wt)

    return(out[order(out$Pval)])
}

run_trials_general <- function(trials, ns, p_ct, fc_ct, norgs, IoD=NULL, IoD_contr, sim_type, modality, cells_per_samp=10000, DA_meth, norm.method='TMM', pca=NULL, io=NULL) {
    trials_out <- lapply(1:trials, function(t) {
        sim <- do_experiment_general(ns=ns, p_ct=p_ct, fc_ct, norgs=norgs, IoD=IoD, IoD_contr=IoD_contr,  sim_type=sim_type, modality=modality, cells_per_samp=10000)
        sim$trial <- t
        if (DA_meth == 'poisson') {
            out <- do_poisson_DA(sim)
        } else if (DA_meth == 'NB') {
            out <- do_NB_DA(sim, norm.method=norm.method)
        } else if (DA_meth == 'Deconv') {
            if (is.null(pca)) stop('pca not provided')
            out <- do_Deconv_DA(sim, io=io, pca=pca)
        } else {
            stop('DA_meth must be either poisson or NB')
        }
        out$trial <- t
        if ((sim_type=="gampois")) out$IoD <- IoD else out$IoD <- NA
        if (!(is.null(IoD_contr))) out$IoD_contr <- IoD_contr else out$IoD_contr <- NA
        out$sim_type <- sim_type
        out$modality <- modality
        out$norgs <- norgs
        out$ns <- ns
        out$nct <- length(p_ct)
        out$cells_per_samp <- cells_per_samp
        return(list(sim=sim, out=out))
    })
    IoD_out <- lapply(trials_out, function(out) out$sim) %>% rbindlist
    IoD_out <- IoD_out[ ,list(mean=mean(ncells), var=var(ncells)), by=c("celltype", "sample", "condition", "IoD", "IoD_contr", "sim_type", "norgs", "ns", "cells_per_samp", "nct", "modality")]
    trials_out <- lapply(trials_out, function(out) out$out) %>% rbindlist
    return(list(trials_out=trials_out, IoD_out=IoD_out))
}

run_trials_gastruloid <- function(to.sample, fc_ct, trials, ns, ng, IoD_contr, IoD=NULL, sim_type, cells_per_samp=10000, DA_meth, norm.method='TMM', pca=NULL, io=NULL) {
    trials_out <- lapply(1:trials, function(t) {
        sim <- do_experiment_gastruloid(to.sample=to.sample, fc_ct, ns=ns, ng=ng, IoD_contr=IoD_contr, IoD=IoD, sim_type=sim_type, cells_per_samp=10000)
        sim$trial <- t
        if (DA_meth == 'poisson') {
            out <- do_poisson_DA(sim)
        } else if (DA_meth == 'NB') {
            out <- do_NB_DA(sim, norm.method=norm.method)
        } else if (DA_meth == 'Deconv') {
            if (is.null(pca)) stop('pca not provided')
            out <- do_Deconv_DA(sim, io=io, pca=pca)
        } else {
            stop('DA_meth must be either poisson or NB')
        }
        out$trial <- t
        if (sim_type=="gampois") out$IoD <- IoD else out$IoD <- NA
        if (!(is.null(IoD_contr))) out$IoD_contr <- IoD_contr else out$IoD_contr <- NA
        out$sim_type <- sim_type
        out$ng <- ng
        out$ns <- ns
        out$cells_per_samp <- cells_per_samp
        return(list(sim=sim, out=out))
    })
    IoD_out <- lapply(trials_out, function(out) out$sim) %>% rbindlist
    IoD_out <- IoD_out[ ,list(mean=mean(ncells), var=var(ncells)), by=c("celltype", "sample", "condition", "IoD", "IoD_contr", "sim_type", "ng", "ns", "cells_per_samp")]
    trials_out <- lapply(trials_out, function(out) out$out) %>% rbindlist
    return(list(trials_out=trials_out, IoD_out=IoD_out))
}

sim_DA_analysis_general <- function(cells_per_samp=10000, fc_ct_list, IoD_seq, nct_seq, n_samps_percond, norg_seq, trials=100, sim_type, modality, p_ct_meth, IoD_contr_seq, dirichlet_param1=NULL, dirichlet_param2=NULL, DA_meth, norm.method='TMM', pca=NULL, io=NULL) {
    
    if (!(sim_type %in% c("gampois", "poisson"))) {
        stop(paste0("sim_type must be one of: ", paste(c("gampois", "poisson"), collapse = ", ")))
    }
    
    if (sim_type != 'gampois') {
        IoD_seq <- c('foo')
    } else if (length(IoD_seq)==0) {
        stop('If sim_type is gampois, IoD seq must be at least of length 1')
    }

    out <- list()
    IoDs <- list()
    i <- 1
    for (nct in nct_seq) {
                    
        if ((p_ct_meth == "uniform") & (modality=="unimodal")) {
            p_ct <- rep(1/nct,nct)
        } else if ((p_ct_meth == "dirichlet") & (modality=="unimodal")) {
            if (is.null(dirichlet_param1)) stop('dirichlet_param1 not provided')
            p_ct <- rdirichlet(1, rep(dirichlet_param1,nct))
        } else if ((p_ct_meth == "dirichlet") & (modality=="bimodal")) {
            if (is.null(dirichlet_param1)) stop('dirichlet_param1 not provided')
            if (is.null(dirichlet_param2)) stop('dirichlet_param2 not provided')
            p_ct_A <- rdirichlet(1, c(rep(dirichlet_param1,ceiling(nct/2)), rep(dirichlet_param2,floor(nct/2))))
            p_ct_B <- rdirichlet(1, c(rep(dirichlet_param2,ceiling(nct/2)), rep(dirichlet_param1,floor(nct/2))))
        } else {
            stop('Unknown p_ct obtainment method')
        }
        
        for (norgs in norg_seq) {
            for (ns in n_samps_percond) {
                for (IoD in IoD_seq) {
                    for (IoD_contr in IoD_contr_seq) {
                        for (fc_ct in fc_ct_list) {

                            if (sim_type == 'gampois') {
                                tmp <- run_trials_general(trials=trials, fc_ct, ns=ns, p_ct=p_ct, norgs=norgs, IoD=IoD, sim_type=sim_type, modality=modality, IoD_contr=IoD_contr, cells_per_samp=10000, DA_meth=DA_meth, norm.method=norm.method, pca=pca, io=io)
                            } else {
                                tmp <- run_trials_general(trials=trials, fc_ct, ns=ns, p_ct=p_ct, norgs=norgs, IoD=NULL, sim_type=sim_type, modality=modality, IoD_contr=IoD_contr, cells_per_samp=10000, DA_meth=DA_meth, norm.method=norm.method, pca=pca, io=io)
                            }
                            out[[i]] <- tmp$trials_out
                            IoDs[[i]] <- tmp$IoD_out
                            out[[i]]$fc_ct <- list(fc_ct)
                            IoDs[[i]]$fc_ct <- list(fc_ct)
                            i <- i + 1

                        }
                    }
                }
            }
        }
    }
    out <- rbindlist(out)
    out[,sig := (Pval < 0.05)]
    out[,discovery := (FDR < 0.1)]
    out[,false_pos := ((Pval < 0.05) & (FDR < 0.1))]
    IoDs <- rbindlist(IoDs)
    IoDs[,IoD_samp := var/mean]
    
    return(list(out=out, IoDs=IoDs))
    
}

do_pca <- function(df_ref, tp, io) {
    
    pca_out <- run_all_pcas(df_ref,
                            npcs=io$npcs[[tp]][["All"]],
                            #npcs=list("mesodermal"=25, "neural"=25),
                            timepoint_values=io$tps_ref[[tp]],
                            classes=io$gastr_type_subset[[tp]]
                           )

    pca <- list()
    pca$X <- pca_out[[2]]
    pca$ct_ordering <- unique(lapply(lapply(pca_out[[4]], function(l) py_to_r(l)), function(X) colnames(X)))[[1]]
    pca$offset <- pca_out[[1]]
    
    return(pca)
}

sim_DA_analysis_gastruloid <- function(g_meta, tp, fc_ct_list, gastr_type_subset=NULL, cells_per_samp=10000, n_samps_percond, ng_seq, IoD_contr_seq, IoD_seq=NULL, trials=100, sim_type, DA_meth, norm.method='TMM', pca=NULL, io=NULL) {
    
    if (!(sim_type %in% c("full", "multinom", "gampois", "poisson", unique(g_meta[!is.na(gastr_type)]$gastr_type)))) {
        stop(paste0("sim_type must be one of: ", paste(c("full", "multinom", "gampois", "poisson", unique(g_meta[!is.na(gastr_type)]$gastr_type)), collapse = ", ")))
    }
    
    if (!(tp %in% unique(g_meta$timepoint))) {
        stop(paste0("tp must be one of: ", paste(unique(g_meta$timepoint), collapse = ", ")))
    }
    
    if (sim_type != 'gampois') {
        IoD_seq <- c('foo')
    } else if (length(IoD_seq)==0) {
        stop('If sim_type is gampois, IoD seq must be at least of length 1')
    }
    
    to.sample <- g_meta[(timepoint == tp) & (cluster != 'unlabelled') & (MULTI_class != 'Negative'), c("cell", "cluster", "MULTI_class", "gastr_type", "timepoint")]
    if (!(is.null(gastr_type_subset))) {
        to.sample <- to.sample[gastr_type %in% gastr_type_subset]
    }
    if (sim_type %in% unique(g_meta[!is.na(gastr_type)]$gastr_type)) {
        to.sample <- to.sample[gastr_type == sim_type]
        sim_type_new <- 'full'
    } else {
        sim_type_new <- sim_type
    }
    to.sample <- to.sample[, sum := .N, by = MULTI_class][, prop := .N, by = c("MULTI_class", "cluster")][, prop := prop/sum][,cell := NULL] %>% unique
    to.sample <- dcast(to.sample, MULTI_class~cluster, value.var = 'prop', fun.aggregate = sum)
    to.sample.gastr_type <- unlist(lapply(to.sample$MULTI_class, function(x) as.character(unique(g_meta[MULTI_class == x]$gastr_type))))
    to.sample$MULTI_class <- NULL
    
    if (is.null(pca) & (DA_meth == 'Deconv')) pca <- do_pca(df_ref=g_meta, tp=tp, io=io)                               
    
    out <- list()
    IoDs <- list()
    i <- 1
    for (ng in ng_seq) {
        for (ns in n_samps_percond) {
            for (IoD_contr in IoD_contr_seq) {
                for (IoD in IoD_seq) {
                    for (fc_ct in fc_ct_list) {
                        
                        if (sim_type == 'gampois') {
                            tmp <- run_trials_gastruloid(to.sample=to.sample, fc_ct, trials=trials, ns=ns, ng=ng, IoD_contr=IoD_contr, IoD=IoD, sim_type=sim_type_new, cells_per_samp=10000, DA_meth=DA_meth, norm.method=norm.method, pca=pca, io=io)
                        } else {
                            tmp <- run_trials_gastruloid(to.sample=to.sample, fc_ct, trials=trials, ns=ns, ng=ng, IoD_contr=IoD_contr, sim_type=sim_type_new, cells_per_samp=10000, DA_meth=DA_meth, norm.method=norm.method, pca=pca, io=io)
                        }
                        out[[i]] <- tmp$trials_out
                        IoDs[[i]] <- tmp$IoD_out
                        out[[i]]$sim_type <- sim_type
                        out[[i]]$fc_ct <- list(fc_ct)
                        IoDs[[i]]$sim_type <- sim_type
                        IoDs[[i]]$fc_ct <- list(fc_ct)
                        i <- i + 1
                        
                    }
                }
            }
        }
    }
    out <- rbindlist(out)
    IoDs <- rbindlist(IoDs)
    out[,sig := (Pval < 0.05)]
    out[,discovery := (FDR < 0.1)]
    out[,false_pos := ((Pval < 0.05) & (FDR < 0.1))]
    IoDs[,IoD_samp := var/mean]
    
    return(list(out=out, IoDs=IoDs))
    
}