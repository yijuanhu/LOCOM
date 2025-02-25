# options(echo=TRUE) # if you want to see commands in output file
# args=(commandArgs(TRUE))
# print(args) 
# 
# if (length(args)==0) {
#     print("No arguments supplied")
    
    n_confounders <- 0
    n_sam <- 152        # 152, 76
    
    causal_type <- 1
    disp1 <- 0.001      # 0.001, 0.01
    depth_fold <- 10
    
    beta <- 1
    
    n_sim <- 1          # number of simulation replicates
    i_seed <- 1         # seed for simulation
    n_cores <- 4
    
# } else {
#     for(i in 1:length(args)) {
#         eval(parse(text=args[[i]]))
#     }
# }
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

library(LDM)
library(dirmult)
library(vegan)
library(GUniFrac)
library(parallel)
library(permute)
library(multtest)
library(LOCOM)


fdr.target <- 0.2 
filter.thresh <- 0.2

filename <- paste("C", n_confounders, "_cau", causal_type, "_n", n_sam, "_beta", beta, "_seed", i_seed,".txt", sep="")

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
    bias.factor1.log[top_taxa] <- 0
    bias.factor2.log[top_taxa] <- 0
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
    
sim = 1
print(sim)

set.seed(i_seed*1000+sim)

# -----------------------------------
# Simulating data
# -----------------------------------

n0.sam <- ceiling(n_sam * 0.5)
n1.sam <- floor(n_sam * 0.5)

Y <- c(rep(0, n0.sam), rep(1, n1.sam))

if (n_confounders > 0) {
    C <- matrix(NA, nrow=n_sam, ncol=n_confounders)
    for (r in 1:n_confounders) {
        C[,r] <- c(runif(n0.sam, -1, 1), runif(n1.sam, 0, 2)) 
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
if (disp > 1e-8) {
    pi.table.sim <- rdirichlet(n = n_sam, (1-disp)/disp*pi_est) 
}
colnames(pi.table.sim) <- c(1:n.otus)

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
otu.table.sim.pool <- otu.table1.sim + otu.table2.sim


# #---------------------------
# # Use the code for real data
# #---------------------------
# 
otu_tab_common <- otu.table1.sim
shotgun_tab_common <- otu.table2.sim
meta_common <- data.frame(int = rep(1, nrow(otu_tab_common)))

var = Y
var_name = "Y"
var_level = sort(unique(var))
var_nlevel = length(var_level)
table(var)

###################################
# Check overdispersion 
###################################

# otu_tab_common_rrf <- Rarefy(otu_tab_common, depth = min(rowSums(otu_tab_common)))$otu.tab.rff # 5834
# shotgun_tab_common_rrf <- Rarefy(shotgun_tab_common, depth = min(rowSums(shotgun_tab_common)))$otu.tab.rff # 8559
# discordance <- matrix(NA, nrow=nrow(otu_tab_common), ncol=ncol(otu_tab_common))
# discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf > 0] <- 0
# discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf == 0] <- 1
# discordance[otu_tab_common_rrf == 0 & shotgun_tab_common_rrf > 0] <- -1
# discordance.col <- apply(discordance, 2, function(x) c(sum(x==1, na.rm=TRUE), sum(x==0, na.rm=TRUE), sum(x==-1, na.rm=TRUE))/sum(table(x)))
# rowMeans(discordance.col, na.rm=TRUE)
# # (sim: 0.01, 0.01, 0.01, 5000, 100000) 0.2716238 0.3919649 0.3364113
# # (ORIGINS) 0.2511531 0.2904511 0.4583958
# # Below are old results
# # # disp_lib=(0.02, 0, 0, 5K, 100K): 0.0146847 0.6380098 0.3473053
# # # disp_lib=(0.02, 0, 0, 10K, 10K): 0.1759458 0.7191579 0.1048963
# # # disp_lib=(0.01, 0, 0, 5K, 100K): 0.01897993 0.62728231 0.35373776
# # # disp_lib=(0.01, 0, 0, 10K, 10K): 0.1963734 0.6929150 0.1107116
# # # disp_lib=(0, 0, 0, 5K, 100K): 0.00771634 0.63458128 0.35770238
# # # disp_lib=(0, 0, 0, 10K, 10K): 0.2528176 0.5995246 0.1476578
# # # disp_lib=(0.01, 0.01, 0.01, 5K, 100K): 0.2838036 0.4139664 0.3022300  (***)
# # # disp_lib=(0.01, 0.01, 0.01, 10K, 10K): 0.3753438 0.3959640 0.2286922
# # # CRTAM data: 0.2891350 0.2901812 0.4206838
# # # CRTAM data (rarefaction): 

###################################
# Scatter plots
###################################

otu_tab_common_freq <- otu_tab_common/rowSums(otu_tab_common)
shotgun_tab_common_freq <- shotgun_tab_common/rowSums(shotgun_tab_common)
o <- order(colMeans(rbind(otu_tab_common_freq, shotgun_tab_common_freq)), decreasing = T)
genus_name_list <- colnames(otu_tab_common)

pdf(paste("./DATA_sim/plot_scatter_sim_n", n_sam, "_d", disp1, "_1-15_36-50.pdf", sep=""), height=18, width=15)
par(mfrow=c(6,5), pty="s", mar=c(3,3,4,0.1), cex.main=1.5, cex.lab=1.5)
for(i in c(1:15, 36:50)){
    genus_name <- paste("Genus ", genus_name_list[o[i]], sep="")
    rho_value <- signif(cor(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]]),2)
    #plot(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]], main = bquote(atop(.(genus_name), rho==.(rho_value))), xlab = "", ylab = "")
    ylab = ifelse(i%%5==1, "SMS RA", "")
    plot(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]], main = bquote(.(genus_name)*","~rho==.(rho_value)), xlab = "", ylab = ylab)
    abline(coef(lm(shotgun_tab_common_freq[,o[i]] ~ otu_tab_common_freq[,o[i]])))
    abline(a = 0, b = 1, col = "red")
}
dev.off()


# ###################################
# # combining data
# ###################################
# 
# n_sam <- nrow(otu_tab_common)
# ID <- 1:n_sam
# ID_comSam <- c(ID, ID)
# var_comSam <- as.factor(c(var, var))
# meta_common_comSam <- rbind(meta_common, meta_common)
# meta_common_comSam$strata <- c(rep(0, n_sam), rep(1, n_sam))
# 
# tab_comSam <- rbind(otu_tab_common, shotgun_tab_common)
# # tab_comTaxa <- cbind(otu_tab_common, shotgun_tab_common)
# tab_aveFreq <- (otu_tab_common/rowSums(otu_tab_common) + shotgun_tab_common/rowSums(shotgun_tab_common))/2
# tab_comCt <- otu_tab_common + shotgun_tab_common
# 
# ###################################
# # Beta diversity
# ###################################
# 
# dist_method = "bray"
# 
# # p.perm.16S <- permanovaFL(otu_tab_common ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
# # p.perm.shotgun <- permanovaFL(shotgun_tab_common ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
# # p.perm.comSam <- permanovaFL(tab_comSam | strata ~ var_comSam, data=meta_common_comSam, dist.method=dist_method, seed=82955,
# #                              cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free")$p.permanova
# # p.perm.aveFreq <- permanovaFL(tab_aveFreq ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
# # p.perm.comCt <- permanovaFL(tab_comCt ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
# # 
# # (p.com.p <- pcauchy(mean(c(tan( (0.5 - p.perm.16S)*pi), tan( (0.5 - p.perm.shotgun)*pi))), lower.tail=F))
# 
# 
# color_16S = rep("blue", n_sam)
# color_16S[which(var==var_level[2])] = "blue"
# color_shotgun = rep("red", n_sam)
# color_shotgun[which(var==var_level[2])] = "red"
# color_com = rep("black", n_sam)
# color_com[which(var==var_level[2])] = "black"
# pch_16S = rep(1, n_sam)
# pch_16S[which(var==var_level[2])] = 16
# pch_shotgun = rep(1, n_sam) 
# pch_shotgun[which(var==var_level[2])] = 16
# pch_com = rep(1, n_sam) 
# pch_com[which(var==var_level[2])] = 16
# 
# dist_16S <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=otu_tab_common, dist.method=dist_method)$adj.dist 
# dist_16S_eigen <- eigen(dist_16S, symmetric=TRUE)
# dist_shotgun <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=shotgun_tab_common, dist.method=dist_method)$adj.dist 
# dist_shotgun_eigen <- eigen(dist_shotgun, symmetric=TRUE)
# dist_comSam <- adjust.data.by.covariates(formula= ~ 1, data=meta_common_comSam, otu.table=tab_comSam, dist.method=dist_method)$adj.dist 
# dist_comSam_eigen <- eigen(dist_comSam, symmetric=TRUE)
# dist_aveFreq <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=tab_aveFreq, dist.method=dist_method)$adj.dist 
# dist_aveFreq_eigen <- eigen(dist_aveFreq, symmetric=TRUE)
# dist_comCt <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=tab_comCt, dist.method=dist_method)$adj.dist 
# dist_comCt_eigen <- eigen(dist_comCt, symmetric=TRUE)
# 
# # dist_shotgun_eigen$vectors[,2] <- - dist_shotgun_eigen$vectors[,2]
# # dist_shotgun_eigen$vectors[,3] <- - dist_shotgun_eigen$vectors[,3]
# 
# pdf(paste("./DATA_sim/plot_beta_", dist_method, "_sim_n", n_sam, "_d", disp1, ".pdf", sep=""), height=12, width=12)
# par(mfrow=c(3,3), pty="s", mar=c(4,3,3,1))
# 
# xrange = c(-0.18, 0.17)
# yrange = c(-0.18, 0.27)
# gap = 0.01
# 
# # scatter plot of PC1
# plot(dist_16S_eigen$vectors[,1], dist_shotgun_eigen$vectors[,1], main="PC1", xlab="16S", ylab="SMS", pch=16)
# abline(0,1)
# # scatter plot of PC2
# plot(dist_16S_eigen$vectors[,2], dist_shotgun_eigen$vectors[,2], main="PC2", xlab="16S", ylab="SMS", pch=16)
# abline(0,1)
# # scatter plot of PC2
# plot(dist_16S_eigen$vectors[,3], dist_shotgun_eigen$vectors[,3], main="PC3", xlab="16S", ylab="SMS", pch=16)
# abline(0,1)
# 
# # 16S
# plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="(a) 16S data", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
# text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
# ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
# legend("topright", legend=c("Y=0", "Y=1"), pch=c(1,16), col=c("blue","blue"))
# # legend("bottomright", legend=paste("p =", signif(p.perm.16S,3)), bty="n")
# 
# # shotgun
# plot(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], main="(b) SMS data", xlab="PC1", ylab="PC2", col=color_shotgun, pch=pch_shotgun, xlim=xrange, ylim=yrange) # here use minus in y-axis to adjust 
# text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
# ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
# legend("topright", legend=c("Y=0", "Y=1"), pch=c(1,16), col=c("red", "red"))
# # legend("bottomright", legend=paste("p =", signif(p.perm.shotgun,3)), bty="n")
# 
# # overlay
# plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="Overlaying figures (a) and (b)", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
# text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
# ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
# points(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], col=color_shotgun, pch=pch_shotgun) # here use minus in y-axis to adjust 
# text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
# ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
# legend("topright", legend=c("16S, Y=0", "16S, Y=1", "SMS, Y=0", "SMS, Y=1"), 
#        pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
# # legend("bottomright", legend=paste("p =", signif(p.com.p,3)), bty="n")
# 
# # combine by row
# plot(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2], main="Combining samples", xlab="PC1", ylab="PC2", 
#      col=c(color_16S, color_shotgun), pch=c(pch_16S, pch_shotgun), xlim=c(-0.13, 0.11), ylim=c(-0.125, 0.2))
# text(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2]+gap, labels=c(ID, ID), col=c(color_16S, color_shotgun), cex=0.5)
# ordiellipse(ord=dist_comSam_eigen$vectors, groups=factor(var_comSam, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
# legend("topright", legend=c("16S, Y=0", "16S, Y=1", "SMS, Y=0", "SMS, Y=1"), 
#        pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
# # legend("bottomright", legend=paste("p =", signif(p.perm.comSam,3)), bty="n")
# 
# # average frequencies
# plot(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2], main="Averaging relative abundance", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
# text(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
# ordiellipse(ord=dist_aveFreq_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
# legend("topright", legend=c("Y=0", "Y=1"), pch=c(1,16), col=c("black", "black"))
# # legend("bottomright", legend=paste("p =", signif(p.perm.aveFreq,3)), bty="n")
# 
# # combine by pooling
# plot(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2], main="Pooling count data", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
# text(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
# ordiellipse(ord=dist_comCt_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
# legend("topright", legend=c("Y=0", "Y=1"), pch=c(1,16), col=c("black", "black"))
# # legend("bottomright", legend=paste("p =", signif(p.perm.comCt,3)), bty="n")
# 
# dev.off()

# prop.presence1 <- colMeans(otu.table1.sim > 0)
# otus.keep1 <- which(prop.presence1 >= filter.thresh)
# otu.table1.sim.filter <- otu.table1.sim[, otus.keep1]
# causal.otus1.filter <- which(otus.keep1 %in% causal.otus)
# noncausal.otus1.filter <-setdiff(1:length(otus.keep1), causal.otus1.filter)
# 
# prop.presence2 <- colMeans(otu.table2.sim > 0)
# otus.keep2 <- which(prop.presence2 >= filter.thresh)
# otu.table2.sim.filter <- otu.table2.sim[, otus.keep2]
# causal.otus2.filter <- which(otus.keep2 %in% causal.otus)
# noncausal.otus2.filter <-setdiff(1:length(otus.keep2), causal.otus2.filter)

prop.presence.pool <- colMeans(otu.table.sim.pool > 0)
otus.keep.pool <- which(prop.presence.pool >= filter.thresh)
# otus.keep.pool <- sort(unique(c(otus.keep1, otus.keep2))) # make comparable
otu.table.sim.pool.filter <- otu.table.sim.pool[, otus.keep.pool]
causal.otus.pool.filter <- which(otus.keep.pool %in% causal.otus)
noncausal.otus.pool.filter <- setdiff(1:length(otus.keep.pool), causal.otus.pool.filter)

bias.factor1.log.pool.filter <- bias.factor1.log[otus.keep.pool]
bias.factor2.log.pool.filter <- bias.factor2.log[otus.keep.pool]

mat1 <- 1:length(otus.keep.pool)
mat2 <- 1:length(otus.keep.pool)
res.locom.new <- locom_com_omni(taxa.table1 = otu.table1.sim[,otus.keep.pool], 
                                taxa.table2 = otu.table2.sim[,otus.keep.pool], 
                                Y1 = Y, Y2 = Y, C1 = C, C2 = C, mat1 = mat1, mat2 = mat2, 
                                cluster.id = factor(rep(1:n_sam, 2)), perm.within.type = "none", perm.between.type = "free",
                                n.perm.max = 50000, fdr.nominal = fdr.target, n.cores = n_cores, n.rej.stop=n_rej_stop)
# n=150, rep=29001

save(res.locom.new, bias.factor1.log.pool.filter, bias.factor2.log.pool.filter, file="sim_n150_d0.001.RData")
# load("sim_n150_d0.001.RData")


dim(res.locom.new$para.est.wc) # 3 212
color <- rep("black", length(bias.factor1.log.pool.filter))
color[bias.factor1.log.pool.filter==-5] <- "red"
color[bias.factor2.log.pool.filter==-5] <- "blue"
color[bias.factor1.log.pool.filter==-5 & bias.factor2.log.pool.filter==-5] <- "purple"


pdf(paste("./DATA_sim/plot_intercept_sim_n", n_sam, "_d", disp1, ".pdf", sep=""), height=4, width=8)
par(mfrow=c(1,2), pty="s", mar=c(5,4,4,2))
plot(res.locom.new$para.est.wc[2,], res.locom.new$para.est.wc[3,], main="New-wc", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16, col=color)
abline(0,1,col="red")
plot(res.locom.new$para.est.we[2,], res.locom.new$para.est.we[3,], main="New-we", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16, col=color)
abline(0,1,col="red")
dev.off()


