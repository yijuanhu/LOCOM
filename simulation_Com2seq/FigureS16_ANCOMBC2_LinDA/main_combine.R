# The following 3 lines are for batch submission in a computing cluster
options(echo=TRUE)
args=(commandArgs(TRUE))
print(args)

#--------------
# Parameters
#--------------

if (length(args)==0) {
  print("No arguments supplied")
  
  # parameters with varying values
  
  n_confounders <- 0    # number of coufounders in Z, taking values 0 or 1
  n_sam <- 100  # number of samples
  
  causal_type <- 1     # M1-M3
  disp1 <- 0.01         # overdispersion para tau1
  depth_fold <- 10
  
  beta <- 2 # exp(beta) in the X-axis of Figure 4
  
  n_sim <- 1   # number of simulation replicates
  i_seed <- 1         # seed for simulation
  n_cores <- 1 # number of cores for parallel computing
  iter_ind <- 50
} else {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

# parameters with fixed values
n_rej_stop <- 100       # a criterion to stop the sequential permutation procedure # the No. of significant results?

bias_sd1 <- 0.5         # sd for simulating experimental bias in dataset 1
bias_sd2 <- 0.5         # sd for simulating experimental bias in dataset 2

disp <- 0.01            # overdispersion para tau
disp1 <- disp1          # overdispersion para tau1
disp2 <- disp1          # overdispersion para tau2
depth.mu1 <- 10000
depth.mu2 <- 10000 * depth_fold

have_bias <- 1
have_diff_bias <- 1
depth.sd1 <- depth.mu1/3
depth.sd2 <- depth.mu2/3
depth.lower <- 2000

library(dirmult) # rdirichlet
library(vegan)
library(abind)   # abind
library(matrixStats) # colMedians
library(parallel)
library(permute)
library(multtest)
library(Rcpp)
library(RcppArmadillo)
library(MicrobiomeStat) # linDA
library(ANCOMBC) # ANCOM-BC2
library(dplyr)
library(TreeSummarizedExperiment)
library(mia)
library(scater)

source("Com2seq_fun_ind.R")
source("LOCOM_fun.R")
sourceCpp("Newton_v2.cpp")

summarize_otu_results <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.2) {

    otu.detected = colnames(qvalue)[which(qvalue < fdr.target)]
    n.otu = length(otu.detected)

    if (n.otu > 0) {
        sen = sum(otu.detected %in% paste0("otu", causal.otus))/length(causal.otus)
        sep = 1 - sum(otu.detected %in% paste0("otu", noncausal.otus))/length(noncausal.otus)
        fdr = n.otu - sum(otu.detected %in% paste0("otu", causal.otus))
        fdr = fdr/n.otu
    } else {
        sen = 0
        sep = 1
        fdr = 0
    }

    out = list(n.otu=n.otu, sen=sen, sep=sep, fdr=fdr)
    return(out)
}


fdr.target <- 0.2        # nominal FDR level
filter.thresh <- 0.2     # filter out rare taxa that have non-zero counts in less than 20% of samples

# File to save simulation results
# current_date <- format(Sys.Date(), format = "%m%d")
filename.our.ind <- paste("Ind","_C", n_confounders, "_cau", causal_type, "_disp", disp1, "_depthFold", depth_fold, "_n", n_sam, "_beta", beta, "_seed", i_seed,".txt", sep="")
filename.linda <- paste("linda","_C", n_confounders, "_cau", causal_type, "_disp", disp1, "_depthFold", depth_fold, "_n", n_sam, "_beta", beta, "_seed", i_seed,".txt", sep="")
filename.ancombc2 <- paste("AncomBC2","_C", n_confounders, "_cau", causal_type, "_disp", disp1, "_depthFold", depth_fold, "_n", n_sam, "_beta", beta, "_seed", i_seed,".txt", sep="")

# Reading relative abundances 'pi' estimated from the throat data
pi_est <- read.table("input_throat/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
n.otus <- length(pi_est) # 856

#------------------------
# causal taxa (associated with T)
#------------------------

if (causal_type == 1) {           # M1
    causal.otus <- c(which(pi_est >= 0.005))[1:20]
} else if (causal_type == 2) {    # M2
    causal.otus <- order(pi_est, decreasing = TRUE)[1:5]
} else if (causal_type == 3) {    # M3
    causal.otus <- sample(which(pi_est >= 0.0005 & pi_est <= 0.001), 50)
}
noncausal.otus <- setdiff(1:n.otus, causal.otus)

beta.otu <- rep(beta, length(causal.otus))
beta.otu.log <- log(beta.otu)

#------------------------
# taxa associated with C
#------------------------

if (n_confounders > 0) {

    betaC <- 1.5 # confounding effect size

    confounding.otus <- matrix(NA, nrow=n_confounders, ncol=5)
    w <- which(pi_est >= 0.005)
    for (r in 1:n_confounders) {
        confounding.otus[r,] <- sort(sample(w, 5))
    }
}

#------------------------
# bias factor
#------------------------

if (have_bias) {

    set.seed(0)
    bias.factor1.log <- rnorm(n.otus, 0, bias_sd1)
    bias.factor2.log <- rnorm(n.otus, 0, bias_sd1)

    if (have_diff_bias) {
        if (causal_type==1) {
            sub1 <- 1:5
            sub2 <- 11:15
        } else if (causal_type==2) {
            sub1 <- c(2, 5)
            sub2 <- c(3, 4)
        } else if (causal_type==3) {
            sub1 <- 1:5
            sub2 <- 11:15
        }

        top5 <- order(pi_est, decreasing = TRUE)[1:5]
        sub1 <- c(causal.otus[sub1], setdiff(sample(noncausal.otus, length(noncausal.otus)/5), top5))
        sub2 <- c(causal.otus[sub2], setdiff(sample(noncausal.otus, length(noncausal.otus)/5), top5))

        bias.factor1.log[sub1] <- -5
        bias.factor2.log[sub2] <- -5
    }
    bias.factor1 <- exp(bias.factor1.log)
    bias.factor2 <- exp(bias.factor2.log)
    top_taxa <- order(pi_est, decreasing = TRUE)[1]
    bias.factor1[top_taxa] <- 1
    bias.factor2[top_taxa] <- 1
}

#-------------------------
# simulation
#-------------------------
for (sim in c(1:n_sim)) {
    print(sim)
    tryCatch({

    set.seed(i_seed*1000+sim)

    # T and C

    n0.sam <- n_sam * 0.5
    n1.sam <- n_sam * 0.5

    Trait <- c(rep(0, n0.sam), rep(1, n1.sam))

    if (n_confounders > 0) {
        C <- matrix(NA, nrow=n_sam, ncol=n_confounders)
        for (r in 1:n_confounders) {
            C[,r] <- c(runif(n0.sam, -1, 1), runif(n1.sam, 0, 2))
        }
    } else {
        C <- NULL
    }

    # Underlying relative abundances
    pi.table.sim <- matrix(rep(pi_est, n_sam), nrow=n_sam, byrow=TRUE)
    if (disp > 1e-8) {
        pi.table.sim <- rdirichlet(n = n_sam, (1-disp)/disp*pi_est)
    }
    colnames(pi.table.sim) <- paste0("otu", c(1:n.otus))
    rownames(pi.table.sim) <- paste0("sam", c(1:n_sam))

    # Introducing the effects of T

    pi.table.sim[, causal.otus] <- pi.table.sim[, causal.otus] * exp(Trait %*% t(beta.otu.log))

    # Introducing the effects of C

    if (n_confounders > 0) {
        for (r in 1:n_confounders) {
            pi.table.sim[, confounding.otus[r,]] <- pi.table.sim[, confounding.otus[r,]] * betaC^C[,r]
        }
    }

    # Introducing experimental bias

    pi.table1.sim <- pi.table.sim
    pi.table2.sim <- pi.table.sim
    if (have_bias == 1) {
        pi.table1.sim <- t(t(pi.table1.sim) * bias.factor1)
        pi.table2.sim <- t(t(pi.table2.sim) * bias.factor2)
    }

    # Normalizing

    pi.table1.sim <- pi.table1.sim/rowSums(pi.table1.sim)
    pi.table2.sim <- pi.table2.sim/rowSums(pi.table2.sim)

    # Generating count data

    otu.table1.sim <- pi.table1.sim # inherit dimension, colnames, rownames
    otu.table2.sim <- pi.table2.sim

    depth1.sim <- rnorm(n_sam, depth.mu1, depth.sd1)
    depth1.sim[depth1.sim < depth.lower] <- depth.lower
    depth1.sim <- round(depth1.sim)

    depth2.sim <- rnorm(n_sam, depth.mu2, depth.sd2)
    depth2.sim[depth2.sim < depth.lower] <- depth.lower
    depth2.sim <- round(depth2.sim)

    for (i in 1:n_sam) {
        if (disp1 < 1e-8) {
            otu.table1.sim[i,] <- rmultinom(1, depth1.sim[i], pi.table1.sim[i,])
        } else {
            otu.table1.sim[i,] <- simPop(J = 1, n = depth1.sim[i], pi = pi.table1.sim[i,], theta = disp1)$data
        }
        if (disp2 < 1e-8) {
            otu.table2.sim[i,] <- rmultinom(1, depth2.sim[i], pi.table2.sim[i,])
        } else {
            otu.table2.sim[i,] <- simPop(J = 1, n = depth2.sim[i], pi = pi.table2.sim[i,], theta = disp2)$data
        }
    }
    
    #################################
    # Com2seq (independent)
    #################################
    
    ####### 16S
    res.locom.16S <- locom(otu.table1.sim, Y=Trait, C = NULL,
                           fdr.nominal = 0.2, seed = NULL,
                           n.perm.max = NULL, n.rej.stop = 100, n.cores = 1,
                           filter.thresh = 0.2, permute = TRUE, verbose = TRUE)

    otu.locom.16S <- summarize_otu_results(res.locom.16S$q.otu,
                                           causal.otus, noncausal.otus)

    ####### Shotgun
    res.locom.shotgun <- locom(otu.table2.sim, Y=Trait, C = NULL,
                               fdr.nominal = 0.2, seed = NULL,
                               n.perm.max = NULL, n.rej.stop = 100, n.cores = 1,
                               filter.thresh = 0.2, permute = TRUE, verbose = TRUE)

    otu.locom.shotgun <- summarize_otu_results(res.locom.shotgun$q.otu,
                                               causal.otus, noncausal.otus)

    ###### com2seq
    res.locom.Com2seq <- Com2seq.ind(table1 = otu.table1.sim,  table2 = otu.table2.sim,
                                     Y1 = Trait, Y2 = Trait, C1 = C, C2 = C, n.perm.max = 50000,
                                     fdr.nominal = fdr.target, n.cores = n_cores, iter_max = iter_ind,
                                     n.rej.stop = n_rej_stop)

    otu.ind.wc <- summarize_otu_results(res.locom.Com2seq$q.taxa.wc,
                                        causal.otus, noncausal.otus)
    tab.com2seq <- c(sim,
                 otu.ind.wc$n.otu,
                 otu.ind.wc$sen,
                 otu.ind.wc$sep,
                 otu.ind.wc$fdr
    )


    # save results
    tab.our.ind = c(t(tab.com2seq), t(unlist(otu.locom.16S)), t(unlist(otu.locom.shotgun)))
    names(tab.our.ind) <- c("sim",
                            "n.otu.com2seq",
                            "sen.com2seq",
                            "sep.com2seq",
                            "fdr.com2seq",
                            "n.otu.16S", "sen.16S","sep.16S", "fdr.16S","n.otu.shotgun", "sen.shotgun","sep.shotgun", "fdr.shotgun")


    if (!file.exists(filename.our.ind)) {
      write.table(t(round(tab.our.ind, 3)), filename.our.ind, append=TRUE, row.names=FALSE,
                  col.names=TRUE, quote=FALSE)
    } else {
      write.table(t(round(tab.our.ind, 3)), filename.our.ind, append=TRUE, row.names=FALSE,
                  col.names=FALSE, quote=FALSE)
    }

    #################################
    # linDA
    #################################
    
    ####### 16S
    # otu.table1.sim.linda = otu.table1.sim
    res.linda.16S <- linda(feature.dat = t(otu.table1.sim), meta.dat = data.frame(Y=Trait), formula = '~Y',
                           alpha = fdr.target,
                           prev.filter = 0.20,
                           feature.dat.type = 'count')
    res.linda.16S.qvalue = t(as.matrix(res.linda.16S$output$Y$padj))
    colnames(res.linda.16S.qvalue) = rownames(res.linda.16S$feature.dat.use)

    otu.linda.16S <- summarize_otu_results(res.linda.16S.qvalue,
                                           causal.otus, noncausal.otus)

    ####### Shotgun

    res.linda.Shotgun <- linda(feature.dat = t(otu.table2.sim), meta.dat = data.frame(Y=Trait), formula = '~Y',
                               alpha = fdr.target,
                               prev.filter = 0.20,
                               feature.dat.type = 'count')
    res.linda.Shotgun.qvalue = t(as.matrix(res.linda.Shotgun$output$Y$padj))
    colnames(res.linda.Shotgun.qvalue) = rownames(res.linda.Shotgun$feature.dat.use)

    otu.linda.Shotgun <- summarize_otu_results(res.linda.Shotgun.qvalue,
                                               causal.otus, noncausal.otus)

    # save results
    tab.linda = c(t(unlist(otu.linda.16S)), t(unlist(otu.linda.Shotgun)))
    names(tab.linda) <- c("n.otu.16S", "sen.16S","sep.16S", "fdr.16S","n.otu.shotgun", "sen.shotgun","sep.shotgun", "fdr.shotgun")


    if (!file.exists(filename.linda)) {
      write.table(t(round(tab.linda, 3)), filename.linda, append=TRUE, row.names=FALSE,
                  col.names=TRUE, quote=FALSE)
    } else {
      write.table(t(round(tab.linda, 3)), filename.linda, append=TRUE, row.names=FALSE,
                  col.names=FALSE, quote=FALSE)
    }
    # 
    #################################
    # ANCOM-BC2
    #################################
    
    ####### 16S
    # Create the tse object
    assays1 = S4Vectors::SimpleList(counts = t(otu.table1.sim))
    smd = data.frame(sample = paste0("sam", seq_len(n_sam)),
                     Y=Trait)
    # smd = data.frame(Y=Trait)
    smd = S4Vectors::DataFrame(smd)
    tse1 = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays1, colData = smd)
    
    #  (SS Filter: pseudo_sens = TRUE, more efficient to control FDR)
    res.ancombc2.16S <-  ancombc2(data = tse1, assay_name = "counts", tax_level = NULL,
                                  fix_formula = "Y", rand_formula = NULL,
                                  p_adj_method = "holm", # default adjustment
                                  alpha = fdr.target,
                                  pseudo_sens = TRUE,
                                  prv_cut = 0.20) 
    res.ancombc2.16S.qvalue = t(as.matrix(res.ancombc2.16S$res$q_Y))
    colnames(res.ancombc2.16S.qvalue) = res.ancombc2.16S$res$taxon
    
    otu.ancombc2.16S  <- summarize_otu_results(res.ancombc2.16S.qvalue,
                                               causal.otus, noncausal.otus)
    
    
    ####### Shotgun 
    assays2 = S4Vectors::SimpleList(counts = t(otu.table2.sim))
    smd = data.frame(sample = paste0("sam", seq_len(n_sam)),
                     Y=Trait)
    smd = S4Vectors::DataFrame(smd)
    tse2 = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays2, colData = smd)    
    #  (SS Filter: pseudo_sens = TRUE, more efficient to control FDR)
    res.ancombc2.Shotgun <- ancombc2(data = tse2, assay_name = "counts", tax_level = NULL,
                                     fix_formula = "Y", rand_formula = NULL,
                                     p_adj_method = "holm", # default adjustment
                                     alpha = fdr.target,
                                     pseudo_sens = TRUE,
                                     prv_cut = 0.20) 
    
    res.ancombc2.Shotgun.qvalue = t(as.matrix(res.ancombc2.Shotgun$res$q_Y))
    colnames(res.ancombc2.Shotgun.qvalue) = res.ancombc2.Shotgun$res$taxon
    
    otu.ancombc2.shotgun  <- summarize_otu_results(res.ancombc2.Shotgun.qvalue,
                                                   causal.otus, noncausal.otus)
    
    # save results
    tab.ancombc2 = c(t(unlist(otu.ancombc2.16S)), t(unlist(otu.ancombc2.shotgun)))
    names(tab.ancombc2) <- c("n.otu.16S", "sen.16S","sep.16S", "fdr.16S","n.otu.shotgun", "sen.shotgun","sep.shotgun", "fdr.shotgun")
    
    if (!file.exists(filename.ancombc2)) {
      write.table(t(round(tab.ancombc2, 3)), filename.ancombc2, append=TRUE, row.names=FALSE,
                  col.names=TRUE, quote=FALSE)
    } else {
      write.table(t(round(tab.ancombc2, 3)), filename.ancombc2, append=TRUE, row.names=FALSE,
                  col.names=FALSE, quote=FALSE)
    }
   
    
  }, error = function(e) {
      message("An error occurred in simulation ", sim, ": ", e$message)
    })
}
