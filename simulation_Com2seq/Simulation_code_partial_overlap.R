options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args) 

if (length(args)==0) {
    print("No arguments supplied")
    
    n_confounders <- 0
    n_sam <- 100
    
    causal_type <- 2
    disp1 <- 0.01
    depth_fold <- 10
    
    beta <- 1.8
    
    n_sim <- 1          # number of simulation replicates
    i_seed <- 1         # seed for simulation
    n_cores <- 4
    
} else {
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}

##### Edited for "general" scenario
n_sam_common <- 40
n_sam_unique1 <- 40
n_sam_unique2 <- 20
pCase_common <- 15/40
pCase_unique1 <- 30/40
pCase_unique2 <- 5/20

i_sam_1 <- c(1:n_sam_unique1, (1:n_sam_common)+n_sam_unique1)
i_sam_2 <- c((1:n_sam_common)+n_sam_unique1, (1:n_sam_unique2)+n_sam_common+n_sam_unique1)
Y <- rep(0, n_sam)
i_sam_Y1 <- c(1:(n_sam_unique1*pCase_unique1),
              1:(n_sam_common*pCase_common)+n_sam_unique1,
              1:(n_sam_unique2*pCase_unique2)+n_sam_common+n_sam_unique1)
Y[i_sam_Y1] <- 1
####################################

n_rej_stop <- 100

bias_sd1 <- 0.5
bias_sd2 <- 0.5

disp <- 0.01
disp1 <- disp1
disp2 <- disp1
depth.mu1 <- 10000
depth.mu2 <- 10000 * depth_fold
 
have_bias <- 1
have_diff_bias <- 1
depth.sd1 <- depth.mu1/3 
depth.sd2 <- depth.mu2/3 
depth.lower <- 2000  

library(dirmult)
library(vegan)
library(parallel)
library(permute)
library(multtest)
library(LOCOM)


fdr.target <- 0.2 
filter.thresh <- 0.2

filename <- paste("C", n_confounders, "_cau", causal_type, "_disp", disp1, "_depthFold", depth_fold, "_n", n_sam, "_beta", beta, "_seed", i_seed,".txt", sep="")

#------------------------
# read otu table 
#------------------------

pi_est <- read.table("input_throat/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
n.otus <- length(pi_est) # 856

#------------------------
# causal taxa
#------------------------
# confounding taxa
#------------------------

if (causal_type == 1) {
    causal.otus <- c(which(pi_est >= 0.005))[1:20]
} else if (causal_type == 2) {
    causal.otus <- order(pi_est, decreasing = TRUE)[1:5]
} else if (causal_type == 3) {
    causal.otus <- sample(which(pi_est >= 0.0005 & pi_est <= 0.001), 50)
} 
noncausal.otus <- setdiff(1:n.otus, causal.otus)

beta.otu <- rep(beta, length(causal.otus))
# beta.otu <- runif(length(causal.otus), 1/beta, beta)
beta.otu.log <- log(beta.otu)

if (n_confounders > 0) {
    
    betaC <- 1.5 # confounding effect
    
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


summarize_otu_results <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.2) {
    
    otu.detected = colnames(qvalue)[which(qvalue < fdr.target)]
    n.otu = length(otu.detected)
    
    if (n.otu > 0) {
        sen = sum(otu.detected %in% causal.otus)/length(causal.otus)
        sep = 1 - sum(otu.detected %in% noncausal.otus)/length(noncausal.otus)
        fdr = n.otu - sum(otu.detected %in% causal.otus)
        fdr = fdr/n.otu
    } else {
        sen = 0
        sep = 1
        fdr = 0
    }
    
    out = list(n.otu=n.otu, sen=sen, sep=sep, fdr=fdr)
    
    return(out)
}

#-------------------------
# simulation
#-------------------------

for (sim in c(1:n_sim)) {
    
    # sim = 1
    print(sim)
    
    set.seed(i_seed*1000+sim)
    
    # -----------------------------------
    # Simulating data
    # -----------------------------------
    
    if (n_confounders > 0) {
        C <- matrix(NA, nrow=n_sam, ncol=n_confounders)
        for (r in 1:n_confounders) {
            C[i_sam_Y1, r] <- runif(length(i_sam_Y1), 0, 2)
            C[-i_sam_Y1, r] <- runif(n_sam - length(i_sam_Y1), -1, 1)
        }
    } else {
        C <- NULL
    }

    depth1.sim <- rnorm(n_sam, depth.mu1, depth.sd1)
    depth1.sim[depth1.sim < depth.lower] <- depth.lower
    depth1.sim <- round(depth1.sim)
    
    depth2.sim <- rnorm(n_sam, depth.mu2, depth.sd2)
    depth2.sim[depth2.sim < depth.lower] <- depth.lower
    depth2.sim <- round(depth2.sim)
    
    pi.table.sim <- matrix(rep(pi_est, n_sam), nrow=n_sam, byrow=TRUE) 
    colnames(pi.table.sim) <- 1:n.otus
    rownames(pi.table.sim) <- paste("sam", 1:n_sam, sep="")
    if (disp > 1e-8) {
        pi.table.sim[1:n_sam, 1:n.otus] <- rdirichlet(n = n_sam, (1-disp)/disp*pi_est) 
    }
    
    # Introducing the same effects of Y
    
    pi.table.sim[, causal.otus] <- pi.table.sim[, causal.otus] * exp(Y %*% t(beta.otu.log))

    # Introducing the same effects of C
    
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
    
    otu.table1.sim <- pi.table1.sim
    otu.table2.sim <- pi.table2.sim
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
    ##### Edited for "general" scenario
    otu.table1.sim <- otu.table1.sim[i_sam_1, colSums(otu.table1.sim>0)>=1] # >=5: perm-s2 has most inflated T1 0.062
    otu.table2.sim <- otu.table2.sim[i_sam_2, colSums(otu.table2.sim>0)>=1]
    dim(otu.table1.sim) # 80 675 (5:325)
    dim(otu.table2.sim) # 60 700 (5:381)

    # Filtering otus
    
    prop.presence1 <- colMeans(otu.table1.sim > 0)
    otus.keep1 <- which(prop.presence1 >= filter.thresh)
    otu.table1.sim.filter <- otu.table1.sim[, otus.keep1]  # 80 96
    causal.otus1.filter <- which(otus.keep1 %in% causal.otus)  # not needed
    noncausal.otus1.filter <-setdiff(1:length(otus.keep1), causal.otus1.filter)
    
    prop.presence2 <- colMeans(otu.table2.sim > 0)
    otus.keep2 <- which(prop.presence2 >= filter.thresh)
    otu.table2.sim.filter <- otu.table2.sim[, otus.keep2]  # 60 104
    causal.otus2.filter <- which(otus.keep2 %in% causal.otus)
    noncausal.otus2.filter <-setdiff(1:length(otus.keep2), causal.otus2.filter)
    
    #####################
    # 16S, shotgun
    #####################
    
    res.locom1 <- locom(otu.table = otu.table1.sim.filter, Y = Y[i_sam_1], C = C[i_sam_1], n.perm.max = 50000, 
                        fdr.nominal = fdr.target, n.cores = n_cores, n.rej.stop=n_rej_stop)	

    res.locom2 <- locom(otu.table = otu.table2.sim.filter, Y = Y[i_sam_2], C = C[i_sam_2], n.perm.max = 50000, 
                        fdr.nominal = fdr.target, n.cores = n_cores, n.rej.stop=n_rej_stop)	

    #####################
    # Com-ct
    #####################
    
    sam.union <- union(rownames(otu.table1.sim), rownames(otu.table2.sim))
    n.sam.union <- length(sam.union)
    taxa.union <- union(colnames(otu.table1.sim), colnames(otu.table2.sim))
    n.taxa.union <- length(taxa.union)
    
    otu_tab_fill0 <- matrix(0, nrow = n.sam.union, ncol = n.taxa.union)
    shotgun_tab_fill0 <- matrix(0, nrow = n.sam.union, ncol = n.taxa.union)
    otu_tab_fill0[match(rownames(otu.table1.sim), sam.union), match(colnames(otu.table1.sim), taxa.union)] <- otu.table1.sim
    shotgun_tab_fill0[match(rownames(otu.table2.sim), sam.union), match(colnames(otu.table2.sim), taxa.union)] <- otu.table2.sim
    
    pooled_tab <- otu_tab_fill0 + shotgun_tab_fill0
    rownames(pooled_tab) <- sam.union
    colnames(pooled_tab) <- taxa.union
    
    pooled_tab_filter <- pooled_tab[, which(colMeans(pooled_tab > 0) >= filter.thresh)]
    dim(pooled_tab_filter) # 100 114
    
    res.locom.pool <- locom(otu.table = pooled_tab_filter, Y = Y, C = C, n.perm.max = 50000, 
                          fdr.nominal = fdr.target, n.cores = n_cores, n.rej.stop=n_rej_stop)
    
    
    ###########
    # New
    ###########
    
    Y1 <- Y[i_sam_1]
    Y2 <- Y[i_sam_2]
    C1 <- NULL
    C2 <- NULL
    res.locom.new3 <-locom_com_omni_general(table1 = otu.table1.sim, table2 = otu.table2.sim, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, n.perm.max = 50000, 
                                           filter_thresh = filter.thresh, fdr.nominal = fdr.target, n.rej.stop = n_rej_stop, n.cores = n_cores)
    if (abs(beta-1) < 1e-6) {
        res.locom.new2 <-locom_com_omni_general(table1 = otu.table1.sim, table2 = otu.table2.sim, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, n.perm.max = 50000, 
                                            filter_thresh = filter.thresh, fdr.nominal = fdr.target, n.rej.stop = n_rej_stop, n.cores = n_cores, perm.strata=2)
        res.locom.new1 <-locom_com_omni_general(table1 = otu.table1.sim, table2 = otu.table2.sim, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, n.perm.max = 50000, 
                                            filter_thresh = filter.thresh, fdr.nominal = fdr.target, n.rej.stop = n_rej_stop, n.cores = n_cores, perm.strata=1)
    } else {
        res.locom.new2 <- NULL
        res.locom.new2$p.global.wc <- NA
        res.locom.new2$p.global.we <- NA
        res.locom.new2$p.global.omni <- NA
        res.locom.new1 <- NULL
        res.locom.new1$p.global.wc <- NA
        res.locom.new1$p.global.we <- NA
        res.locom.new1$p.global.omni <- NA
    }
    
    ###########
    # p.combine
    ###########
    
    name1 <- colnames(res.locom1$p.otu)
    name2 <- colnames(res.locom2$p.otu)
    common.name <- intersect(name1, name2)
    j.mat1 <- match(common.name, name1)
    j.mat2 <- match(common.name, name2)
    
    p.both <- rbind(res.locom1$p.otu[j.mat1], res.locom2$p.otu[j.mat2])
    p.comp <- pcauchy(apply( tan( (0.5 - p.both)*pi ), 2, mean), lower.tail = F)
    p.comp.HM <- 2/colSums(1/p.both)
    p.comp.name <- common.name
    if (length(res.locom1$p.otu[-j.mat1]) > 0) {
        p.comp <- c(p.comp, res.locom1$p.otu[-j.mat1])
        p.comp.HM <- c(p.comp.HM, res.locom1$p.otu[-j.mat1])
        p.comp.name <- c(p.comp.name, name1[-j.mat1])
    }
    if (length(res.locom2$p.otu[-j.mat2]) > 0) {
        p.comp <- c(p.comp, res.locom2$p.otu[-j.mat2])
        p.comp.HM <- c(p.comp.HM, res.locom2$p.otu[-j.mat2])
        p.comp.name <- c(p.comp.name, name2[-j.mat2])
    }

    p.comp <- matrix(p.comp, nrow=1)
    p.comp.HM <- matrix(p.comp.HM, nrow=1)
    q.comp <- matrix(p.adjust(p.comp, method ="BH"), nrow=1)
    q.comp.HM <- matrix(p.adjust(p.comp.HM, method ="BH"), nrow=1)
    colnames(p.comp) <- p.comp.name
    colnames(q.comp) <- p.comp.name
    colnames(p.comp.HM) <- p.comp.name
    colnames(q.comp.HM) <- p.comp.name
    
    p.global.comp <- pcauchy(mean(tan( (0.5 - p.comp)*pi )), lower.tail = F)
    p.globalHM.comp <- length(p.comp.HM)/sum(1/p.comp.HM)
    
    
    length(res.locom1$q.otu) # 93
    length(res.locom2$q.otu) # 111
    length(q.comp) # 134
    length(res.locom.new3$q.taxa.omni) # 136
    
    # -----------------------------------
    # Summarizing results
    # -----------------------------------
    
    otu.new.wc <- summarize_otu_results(res.locom.new3$q.taxa.wc, causal.otus, noncausal.otus)
    otu.new.we <- summarize_otu_results(res.locom.new3$q.taxa.we, causal.otus, noncausal.otus)
    otu.new.omni <- summarize_otu_results(res.locom.new3$q.taxa.omni, causal.otus, noncausal.otus)
    otu.comp <- summarize_otu_results(q.comp, causal.otus, noncausal.otus)
    otu.comp.HM <- summarize_otu_results(q.comp.HM, causal.otus, noncausal.otus)
    otu.pool <- summarize_otu_results(res.locom.pool$q.otu, causal.otus, noncausal.otus)
    otu1 <- summarize_otu_results(res.locom1$q.otu, causal.otus, noncausal.otus)
    otu2 <- summarize_otu_results(res.locom2$q.otu, causal.otus, noncausal.otus)
    
    tab <- c(res.locom.new3$p.global.wc, res.locom.new3$p.global.we, res.locom.new3$p.global.omni,
             res.locom.new2$p.global.wc, res.locom.new2$p.global.we, res.locom.new2$p.global.omni,
             res.locom.new1$p.global.wc, res.locom.new1$p.global.we, res.locom.new1$p.global.omni,
             p.global.comp, p.globalHM.comp, res.locom.pool$p.global, 
             res.locom1$p.global, res.locom2$p.global, 
             otu.new.wc$n.otu, otu.new.wc$sen, otu.new.wc$fdr,
             otu.new.we$n.otu, otu.new.we$sen, otu.new.we$fdr,
             otu.new.omni$n.otu, otu.new.omni$sen, otu.new.omni$fdr,
             otu.comp$n.otu, otu.comp$sen, otu.comp$fdr, 
             otu.comp.HM$n.otu, otu.comp.HM$sen, otu.comp.HM$fdr, 
             otu.pool$n.otu, otu.pool$sen, otu.pool$fdr, 
             otu1$n.otu, otu1$sen, otu1$fdr,  
             otu2$n.otu, otu2$sen, otu2$fdr)
    
    write.table(t(round(tab, 5)), filename, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    sort(causal.otus[causal.otus %in% colnames(q.comp)])
    # 1  35  59  83  91  98 110 113 120 152 362 390 412 442 562 566 639 696 698 699 742 801 809 842 846
    sort(causal.otus[causal.otus %in% colnames(res.locom.new3$q.taxa.omni)])
    # 1  35  59  83  91  98 110 113 120 152 362 390 412 442 562 566 639 696 698 699 742 801 809 842 846
    all(colnames(q.comp) %in% colnames(res.locom.new3$q.taxa.omni))
    
    length(q.comp)
    length(res.locom.new3$q.taxa.omni)
}


dat <- read.table(filename)
dim(dat)
colMeans(dat)

