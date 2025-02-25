#-----------------
# common taxa
#-----------------

common_taxa <- intersect(colnames(otu_tab_commonSam), colnames(shotgun_tab_commonSam))
length(common_taxa) # 125

# 16S 
otu_tab_common <- otu_tab_commonSam[, match(common_taxa, colnames(otu_tab_commonSam))]
dim(otu_tab_common) # 152 125
otu_tab_common[1:3, 1:3]
#                     g__Porphyromonas g__Akkermansia g__Coprococcus
# 11808.3069.SALB.W2V1              714              0              0
# 11808.3070.SALB.W2V1               15              0              0
# 11808.3115.SALB.W2V1               60              0              0
summary(rowSums(otu_tab_common))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5834   19328   23376   23450   27047   46124 

# shotgun
shotgun_tab_common <- shotgun_tab_commonSam[, match(common_taxa, colnames(shotgun_tab_commonSam))]
dim(shotgun_tab_common) # 152 125
shotgun_tab_common[1:3, 1:3]
#                     g__Porphyromonas g__Akkermansia g__Coprococcus
# 11808.3069.SALB.W2V1             1500              2              0
# 11808.3070.SALB.W2V1              646              0              0
# 11808.3115.SALB.W2V1               24              0              1
summary(rowSums(shotgun_tab_common))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8559   57003  137656  172800  257229 1100108 

###################################
# var
###################################

var = as.factor(meta_common$prediabetes)
var_name = "Prediabetes"  
var_level = sort(unique(var))
var_nlevel = length(var_level)
table(var)
# var
# Not Prediabetes 
# 105          47 

###################################
# Check overdispersion 
###################################
 
otu_tab_common_rrf <- Rarefy(otu_tab_common, depth = min(rowSums(otu_tab_common)))$otu.tab.rff # 5834
shotgun_tab_common_rrf <- Rarefy(shotgun_tab_common, depth = min(rowSums(shotgun_tab_common)))$otu.tab.rff # 8559
discordance <- matrix(NA, nrow=nrow(otu_tab_common), ncol=ncol(otu_tab_common))
discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf > 0] <- 0
discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf == 0] <- 1
discordance[otu_tab_common_rrf == 0 & shotgun_tab_common_rrf > 0] <- -1
discordance.col <- apply(discordance, 2, function(x) c(sum(x==1, na.rm=TRUE), sum(x==0, na.rm=TRUE), sum(x==-1, na.rm=TRUE))/sum(table(x)))
rowMeans(discordance.col, na.rm=TRUE)
#  0.2454932 0.2856715 0.4688353


###################################
# Scatter plots 
###################################

otu_tab_common_freq <- otu_tab_common/rowSums(otu_tab_common)
shotgun_tab_common_freq <- shotgun_tab_common/rowSums(shotgun_tab_common)
o <- order(colMeans(rbind(otu_tab_common_freq, shotgun_tab_common_freq)), decreasing = T)
genus_name_list <- colnames(otu_tab_common)

pdf("plot_scatter_Prediabetes_1-15_36-50.pdf", height=18, width=15)
par(mfrow=c(6,5), pty="s", mar=c(3,3,4,0.1), cex.main=1.5, cex.lab=1.5)
for(i in c(1:15, 36:50)){
  genus_name <- unlist(strsplit(genus_name_list[o[i]], "_"))[3]
  rho_value <- signif(cor(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]]),2)
  #plot(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]], main = bquote(atop(.(genus_name), rho==.(rho_value))), xlab = "", ylab = "")
  ylab = ifelse (i %% 5 == 1, "SMS RA", "")
  plot(otu_tab_common_freq[,o[i]], shotgun_tab_common_freq[,o[i]], main = bquote(.(genus_name)*","~rho==.(rho_value)), xlab = "", ylab = ylab)
  abline(coef(lm(shotgun_tab_common_freq[,o[i]] ~ otu_tab_common_freq[,o[i]])))
  abline(a = 0, b = 1, col = "red")
}
dev.off()


###################################
# combining data
###################################

n_sam <- nrow(otu_tab_common)
ID <- 1:n_sam
ID_comSam <- c(ID, ID)
var_comSam <- as.factor(c(var, var))
meta_common_comSam <- rbind(meta_common, meta_common)
meta_common_comSam$strata <- c(rep(0, n_sam), rep(1, n_sam))

tab_comSam <- rbind(otu_tab_common, shotgun_tab_common)
# tab_comTaxa <- cbind(otu_tab_common, shotgun_tab_common)
tab_aveFreq <- (otu_tab_common/rowSums(otu_tab_common) + shotgun_tab_common/rowSums(shotgun_tab_common))/2
tab_comCt <- otu_tab_common + shotgun_tab_common

###################################
# Alpha diversity
###################################

alpha_measure <- "Shannon"  # Shannon, log Chao1

# 16S
alpha_16S <- if (alpha_measure=="Shannon") diversity(otu_tab_common, index = "shannon", MARGIN = 1) else log(estimateR(otu_tab_common)[2,])
alpha_16S_r <- lm(alpha_16S~0, data=meta_common)$residuals 
# shotgun
alpha_shotgun <- if (alpha_measure=="Shannon") diversity(shotgun_tab_common, index = "shannon", MARGIN = 1) else log(estimateR(shotgun_tab_common)[2,])
alpha_shotgun_r <- lm(alpha_shotgun~0, data=meta_common)$residuals 
# comSam
alpha_comSam <- if (alpha_measure=="Shannon") diversity(tab_comSam, index = "shannon", MARGIN = 1) else log(estimateR(tab_comSam)[2,])
alpha_comSam_r <- lm(alpha_comSam~0, data=meta_common_comSam)$residuals 
alpha_comSam_r <- matrix(alpha_comSam_r, ncol=1)
pvalue.comSam <- signif(ldm(alpha_comSam_r | strata ~ var_comSam, data=meta_common_comSam, seed=82955, dist.method="euclidean",
                            cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free",
                            scale.otu.table=FALSE, freq.scale.only=TRUE)$p.global.freq, 2)
# aveFreq
alpha_aveFreq <- if (alpha_measure=="Shannon") diversity(tab_aveFreq, index = "shannon", MARGIN = 1) else log(estimateR(tab_aveFreq)[2,])
alpha_aveFreq_r <- lm(alpha_aveFreq~0, data=meta_common)$residuals 
# comCt
alpha_comCt <- if (alpha_measure=="Shannon") diversity(tab_comCt, index = "shannon", MARGIN = 1) else log(estimateR(tab_comCt)[2,])
alpha_comCt_r <- lm(alpha_comCt~0, data=meta_common)$residuals 


pdf(paste("plot_alpha_", alpha_measure, "_", var_name, ".pdf", sep=""), height=5.5, width=8)
par(mfrow=c(2,3), pty="s", mar=c(4.5,3,4,0.1))

# scatter plot
plot(alpha_16S_r, alpha_shotgun_r, main = alpha_measure, xlab = "16S data", ylab = "SMS data")
abline(a = 0, b = 1, col = "red")
# 16S
boxplot(alpha_16S_r ~ var, main="16S data", xlab="", ylab=alpha_measure) 
pvalue.16 = signif(wilcox.test(alpha_16S_r ~ var)$p.value, 2)
legend(x="bottomright", legend=paste("p =", pvalue.16), bty="n")
# shotgun
boxplot(alpha_shotgun_r ~ var, main="SMS data", xlab="", ylab=alpha_measure) 
pvalue.shotgun = signif(wilcox.test(alpha_shotgun_r ~ var)$p.value, 2)
legend(x="bottomright", legend=paste("p =", pvalue.shotgun), bty="n")
# com-p
(p.cauchy <- signif(pcauchy(mean(c(tan( (0.5 - pvalue.16)*pi), tan( (0.5 - pvalue.shotgun)*pi))), lower.tail=F), 2))
# 0.028
# comSam
boxplot(alpha_comSam_r ~ var_comSam, main="Combined observations", xlab="", ylab=alpha_measure) 
legend(x="bottomright", legend=paste("p =", pvalue.comSam), bty="n")
# aveFreq
boxplot(alpha_aveFreq_r ~ var, main="Averaged relative abundance", xlab="", ylab=alpha_measure) 
pvalue = signif(wilcox.test(alpha_aveFreq_r ~ var)$p.value, 2)
legend(x="bottomright", legend=paste("p =", pvalue), bty="n")
# comCt
boxplot(alpha_comCt_r ~ var, main="Pooled count data", xlab="", ylab=alpha_measure) 
pvalue = signif(wilcox.test(alpha_comCt_r ~ var)$p.value, 2)
legend(x="bottomright", legend=paste("p =", pvalue), bty="n")

dev.off()


###################################
# Beta diversity
###################################

dist_method = "bray"

p.perm.16S <- permanovaFL(otu_tab_common ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.shotgun <- permanovaFL(shotgun_tab_common ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.comSam <- permanovaFL(tab_comSam | strata ~ var_comSam, data=meta_common_comSam, dist.method=dist_method, seed=82955,
                             cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free")$p.permanova
p.perm.aveFreq <- permanovaFL(tab_aveFreq ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.comCt <- permanovaFL(tab_comCt ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova

(p.com.p <- pcauchy(mean(c(tan( (0.5 - p.perm.16S)*pi), tan( (0.5 - p.perm.shotgun)*pi))), lower.tail=F))
# 0.1638216
pcauchy(mean(c(tan( (0.5 - p.perm.comSam)*pi), tan( (0.5 - p.perm.aveFreq)*pi))), lower.tail=F)
# 0.2645166

color_16S = rep("blue", n_sam)
color_16S[which(var==var_level[2])] = "blue"
color_shotgun = rep("red", n_sam)
color_shotgun[which(var==var_level[2])] = "red"
color_com = rep("black", n_sam)
color_com[which(var==var_level[2])] = "black"
pch_16S = rep(1, n_sam)
pch_16S[which(var==var_level[2])] = 16
pch_shotgun = rep(1, n_sam) 
pch_shotgun[which(var==var_level[2])] = 16
pch_com = rep(1, n_sam) 
pch_com[which(var==var_level[2])] = 16

dist_16S <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=otu_tab_common, dist.method=dist_method)$adj.dist 
dist_16S_eigen <- eigen(dist_16S, symmetric=TRUE)
dist_shotgun <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=shotgun_tab_common, dist.method=dist_method)$adj.dist 
dist_shotgun_eigen <- eigen(dist_shotgun, symmetric=TRUE)
dist_comSam <- adjust.data.by.covariates(formula= ~ 1, data=meta_common_comSam, otu.table=tab_comSam, dist.method=dist_method)$adj.dist 
dist_comSam_eigen <- eigen(dist_comSam, symmetric=TRUE)
dist_aveFreq <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=tab_aveFreq, dist.method=dist_method)$adj.dist 
dist_aveFreq_eigen <- eigen(dist_aveFreq, symmetric=TRUE)
dist_comCt <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=tab_comCt, dist.method=dist_method)$adj.dist 
dist_comCt_eigen <- eigen(dist_comCt, symmetric=TRUE)

dist_shotgun_eigen$vectors[,2] <- - dist_shotgun_eigen$vectors[,2]
dist_shotgun_eigen$vectors[,3] <- - dist_shotgun_eigen$vectors[,3]
dist_comSam_eigen$vectors[,2] <- - dist_comSam_eigen$vectors[,2]
dist_aveFreq_eigen$vectors[,2] <- - dist_aveFreq_eigen$vectors[,2]
dist_comCt_eigen$vectors[,2] <- - dist_comCt_eigen$vectors[,2]

pdf(paste("plot_beta_", dist_method, "_", var_name, ".pdf", sep=""), height=12, width=12)
par(mfrow=c(3,3), pty="s", mar=c(4,3,3,1))

xrange = c(-0.18, 0.17)
yrange = c(-0.18, 0.27)
gap = 0.01

# scatter plot of PC1
plot(dist_16S_eigen$vectors[,1], dist_shotgun_eigen$vectors[,1], main="PC1", xlab="16S", ylab="SMS", pch=16)
abline(0,1)
# scatter plot of PC2
plot(dist_16S_eigen$vectors[,2], dist_shotgun_eigen$vectors[,2], main="PC2", xlab="16S", ylab="SMS", pch=16)
abline(0,1)
# scatter plot of PC2
plot(dist_16S_eigen$vectors[,3], dist_shotgun_eigen$vectors[,3], main="PC3", xlab="16S", ylab="SMS", pch=16)
abline(0,1)

# 16S
plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="(a) 16S data", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
legend("topright", legend=c("Not", "Prediabetes"), pch=c(1,16), col=c("blue","blue"))
legend("bottomright", legend=paste("p =", signif(p.perm.16S,3)), bty="n")

# shotgun
plot(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], main="(b) SMS data", xlab="PC1", ylab="PC2", col=color_shotgun, pch=pch_shotgun, xlim=xrange, ylim=yrange) # here use minus in y-axis to adjust 
text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
legend("topright", legend=c("Not", "Prediabetes"), pch=c(1,16), col=c("red", "red"))
legend("bottomright", legend=paste("p =", signif(p.perm.shotgun,3)), bty="n")

# overlay
plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="Overlaying figures (a) and (b)", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
points(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], col=color_shotgun, pch=pch_shotgun) # here use minus in y-axis to adjust 
text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
legend("topright", legend=c("16S, Not", "16S, Prediabetes", "SMS, Not", "SMS, Prediabetes"), 
       pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
legend("bottomright", legend=paste("p =", signif(p.com.p,3)), bty="n")

# combine by row
plot(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2], main="Combining samples", xlab="PC1", ylab="PC2", 
     col=c(color_16S, color_shotgun), pch=c(pch_16S, pch_shotgun), xlim=c(-0.13, 0.11), ylim=c(-0.125, 0.2))
text(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2]+gap, labels=c(ID, ID), col=c(color_16S, color_shotgun), cex=0.5)
ordiellipse(ord=dist_comSam_eigen$vectors, groups=factor(var_comSam, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("topright", legend=c("16S, Not", "16S, Prediabetes", "SMS, Not", "SMS, Prediabetes"), 
       pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
legend("bottomright", legend=paste("p =", signif(p.perm.comSam,3)), bty="n")

# average frequencies
plot(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2], main="Averaging relative abundance", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
text(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
ordiellipse(ord=dist_aveFreq_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("topright", legend=c("Not", "Prediabetes"), pch=c(1,16), col=c("black", "black"))
legend("bottomright", legend=paste("p =", signif(p.perm.aveFreq,3)), bty="n")

# combine by pooling
plot(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2], main="Pooling count data", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
text(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
ordiellipse(ord=dist_comCt_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("topright", legend=c("Not", "Prediabetes"), pch=c(1,16), col=c("black", "black"))
legend("bottomright", legend=paste("p =", signif(p.perm.comCt,3)), bty="n")

dev.off()


###################################
# LOCOM
###################################

source("../LOCOM_com_omni_fun_general.R")

otu_tab_common_filter <- otu_tab_common[, which(colMeans(otu_tab_common > 0) >= filter_thresh)] 
dim(otu_tab_common_filter) # 152 54
shotgun_tab_common_filter <- shotgun_tab_common[, which(colMeans(shotgun_tab_common > 0) >= filter_thresh)] 
dim(shotgun_tab_common_filter) # 152 90
var_num <- as.numeric(var)

# M1: no confounder
C <- NULL

fdr_target <- 0.1

#-------------
# 16S
#-------------
res.locom.16S <- locom(otu.table = otu_tab_common_filter, Y = var_num, C = C, ref.otu = 22, # "g__Prevotella" 22
                       fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.16S$p.global # (20%) 0.07639236 (10%) 0.06888506
length(res.locom.16S$p.otu) # 53
res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target] 
# g__Haemophilus g__Campylobacter  g__Butyrivibrio g__Chelonobacter 
# 0.1205879        0.1025897        0.1930307        0.1205879 

# g__Campylobacter 
# 0.09299483 
# res.locom.16S$ref
# g__Haemophilus 
# 19 

#-------------
# shotgun
#-------------
res.locom.shotgun <- locom(otu.table = shotgun_tab_common_filter, Y = var_num, C = C, ref.otu = 34, # "g__Prevotella"
                           fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.shotgun$p.global # (20%) 0.1666833 (10%) 0.1699057
length(res.locom.shotgun$p.otu) # 90
res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target] # numeric(0)
    
#-------------
# Com-p 
#-------------
int.name <- sort(intersect(colnames(res.locom.16S$p.otu), colnames(res.locom.shotgun$p.otu)))
j.mat1 <- match(int.name, colnames(res.locom.16S$p.otu))
j.mat2 <- match(int.name, colnames(res.locom.shotgun$p.otu))

p.both <- rbind(res.locom.16S$p.otu[j.mat1], res.locom.shotgun$p.otu[j.mat2])
p.comp <- matrix(pcauchy(apply(tan((0.5-p.both)*pi), 2, mean), lower.tail = F), nrow = 1)
p.comp.HM <- 2/colSums(1/p.both)
p.comp.name <- int.name
if (length(res.locom.16S$p.otu[-j.mat1])>0) {
    p.comp <- c(p.comp, res.locom.16S$p.otu[-j.mat1])
    p.comp.HM <- c(p.comp.HM, res.locom.16S$p.otu[-j.mat1])
    p.comp.name <-c(p.comp.name, colnames(res.locom.16S$p.otu)[-j.mat1])
}
if (length(res.locom.shotgun$p.otu[-j.mat2])>0) {
    p.comp <- c(p.comp, res.locom.shotgun$p.otu[-j.mat2])
    p.comp.HM <- c(p.comp.HM, res.locom.shotgun$p.otu[-j.mat2])
    p.comp.name <-c(p.comp.name, colnames(res.locom.shotgun$p.otu)[-j.mat2])
}
# taxa
p.comp <- matrix(p.comp, nrow=1)
p.comp.HM <- matrix(p.comp.HM, nrow=1)
q.comp <- matrix(p.adjust(p.comp, method ="BH"), nrow = 1)
q.comp.HM <- matrix(p.adjust(p.comp.HM, method ="BH"), nrow = 1)
colnames(p.comp) <- p.comp.name
colnames(q.comp) <- p.comp.name
colnames(q.comp.HM) <- p.comp.name
colnames(p.comp.HM) <- p.comp.name

q.comp[,which(q.comp < fdr_target)]  # numeric(0)
q.comp.HM[,which(q.comp.HM < fdr_target)]  # numeric(0)

# global p-value
pcauchy(mean(tan( (0.5-p.comp)*pi) ), lower.tail = F) # 0.110133  0.106142
length(p.comp.HM)/sum(1/p.comp.HM) # 0.07027911   0.06841047

length(q.comp) # 96

#-------------
# Com-ct
#-------------
tab_comCt <- otu_tab_common + shotgun_tab_common
otus.keep.count <- which(colMeans(tab_comCt > 0) >= filter_thresh) 
length(otus.keep.count) # 98
tab_comCt_filter <- tab_comCt[, otus.keep.count]

res.locom.ct <- locom(otu.table = tab_comCt_filter, Y = var_num, C = C,
                      fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.ct$p.global # 0.02259774  0.02259774
res.locom.ct$q.otu[,res.locom.ct$q.otu < fdr_target] 
# g__Campylobacter  g__Butyrivibrio       g__Gemella 
# 0.02939706       0.02449755       0.02449755 

length(res.locom.ct$p.otu)

#-------------
# Com-new
#-------------
length(otus.keep.count) # 98
mat1 <- 1:length(otus.keep.count)
mat2 <- 1:length(otus.keep.count)
res.locom.new <- Com2seq(taxa.table1 = otu_tab_common[,otus.keep.count], taxa.table2 = shotgun_tab_common[,otus.keep.count], Y1 = var_num, Y2 = var_num, C1 = C, C2 = C, mat1 = mat1, mat2 = mat2, ref.taxa = 35, # 35
                              cluster.id=factor(rep(1:n_sam, 2)), perm.within.type="none", perm.between.type="free",
                              fdr.nominal = fdr_target, seed=82955, n.perm.max=50000) # permutation: 18001
res.locom.new$p.global.wc # 0.01979802  0.02605126
res.locom.new$q.taxa.wc[,res.locom.new$q.taxa.wc < fdr_target] 
# g__Pseudoalteromonas        g__Mannheimia           g__Kocuria 
# 0.14138586           0.14138586           0.04246242 
# g__Campylobacter        g__Parvimonas       g__Pasteurella 
# 0.02449755           0.19271406           0.19230577 
# g__Butyrivibrio           g__Gemella       g__Lonepinella 
# 0.14138586           0.02449755           0.12983702 

# g__Kocuria g__Campylobacter       g__Gemella 
# 0.06361069       0.03352455       0.03352455 
length(res.locom.new$q.taxa.wc[,res.locom.new$q.taxa.wc < fdr_target] ) # 9  3

res.locom.new$p.global.we # 0.04979502  0.04357665
res.locom.new$q.taxa.we[,res.locom.new$q.taxa.we < fdr_target] 
# g__Pseudoalteromonas       g__Actinomyces        g__Mannheimia 
# 0.07512582           0.07512582           0.16266373 
# g__Capnocytophaga       g__Lonepinella 
# 0.07512582           0.08084192 

# g__Pseudoalteromonas       g__Actinomyces    g__Capnocytophaga       g__Lonepinella 
# 0.06704910           0.06704910           0.06704910           0.07607494 
length(res.locom.new$q.taxa.we[,res.locom.new$q.taxa.we < fdr_target] ) # 5  4

res.locom.new$p.global.omni # 0.02689731  0.0345245
res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < fdr_target] 
# g__Pseudoalteromonas       g__Actinomyces        g__Mannheimia 
# 0.06532680           0.05552778           0.18944772 
# g__Capnocytophaga           g__Kocuria     g__Campylobacter 
# 0.06532680           0.05879412           0.04899510 
# g__Butyrivibrio           g__Gemella       g__Lonepinella 
# 0.18944772           0.04899510           0.08119188 

# g__Pseudoalteromonas       g__Actinomyces    g__Capnocytophaga           g__Kocuria 
# 0.06361069           0.04469940           0.06361069           0.06361069 
# g__Campylobacter           g__Gemella       g__Lonepinella 
# 0.04469940           0.04469940           0.07736435 
length(res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < fdr_target] ) # 9  7

res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < 0.3] 


save(res.locom.16S, res.locom.shotgun, res.locom.ct, res.locom.new, file="res.overlappingData.FDR10.RData")
load("res.overlappingData.FDR10.RData")

pdf(paste("plot_intercept_", var_name, ".pdf", sep=""), height=4, width=8)
par(mfrow=c(1,2), pty="s", mar=c(5,4,4,2))
plot(res.locom.new$para.est.wc[2,], res.locom.new$para.est.wc[3,], main="New-count", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16, xlim=c(-10, 10))
abline(0,1,col="red")
plot(res.locom.new$para.est.we[2,], res.locom.new$para.est.we[3,], main="New-equal", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16, xlim=c(-10, 10))
abline(0,1,col="red")
dev.off()


#################
# summarizing
#################

taxa_union <- c("g__Gemella", "g__Actinomyces", "g__Campylobacter", "g__Kocuria", "g__Pseudoalteromonas", "g__Capnocytophaga", "g__Lonepinella", "g__Butyrivibrio", "g__Chelonobacter")
taxa_union_16S <- c("g__Actinomyces", "g__Campylobacter", "g__Capnocytophaga", "g__Butyrivibrio", "g__Chelonobacter")
signif(res.locom.new$p.taxa.omni[, taxa_union],3)
signif(res.locom.new$p.taxa.we[, taxa_union],3)
signif(res.locom.new$p.taxa.wc[, taxa_union],3)
signif(res.locom.ct$p.otu[, taxa_union],3)
signif(p.comp[, taxa_union],3)
signif(p.comp.HM[, taxa_union],3)
signif(res.locom.16S$p.otu[, taxa_union_16S],3)
signif(res.locom.shotgun$p.otu[, taxa_union],3)

signif(res.locom.new$q.taxa.omni[, taxa_union],3)
signif(res.locom.new$q.taxa.we[, taxa_union],3)
signif(res.locom.new$q.taxa.wc[, taxa_union],3)
signif(res.locom.ct$q.otu[, taxa_union],3)
signif(q.comp[, taxa_union],3)
signif(q.comp.HM[, taxa_union],3)
signif(res.locom.16S$q.otu[, taxa_union_16S],3)
signif(res.locom.shotgun$q.otu[, taxa_union],3)


# pdf("plot_histogram.pdf", height=8, width=12)
# par(mfrow=c(2,3), pty="s")
# hist(res.locom.16S$p.otu, main="16S data", xlab="p-value")
# hist(res.locom.shotgun$p.otu, main="shotgun data", xlab="p-value")
# hist(p.comp, main="Com-p", xlab="p-value")
# hist(res.locom.new$p.taxa.wc, main="Com-new-wc", xlab="p-value")
# hist(res.locom.new$p.taxa.we, main="Com-new-we", xlab="p-value")
# hist(res.locom.new$p.taxa.omni, main="Com-new-omni", xlab="p-value")
# hist(res.locom.ct$p.otu, main="Com-ct", xlab="p-value")
# dev.off()

# sort(res.locom.16S$p.otu)[1:5]
# sort(res.locom.shotgun$p.otu)[1:5]
# sort(p.comp)[1:5]
# sort(res.locom.new$p.taxa.wc)[1:5]
# sort(res.locom.new$p.taxa.we)[1:5]
# sort(res.locom.new$p.taxa.omni)[1:5]
# sort(res.locom.ct$p.otu)[1:5]


# taxa_union
# [1] "g__Gemella"           "g__Actinomyces"       "g__Campylobacter"     "g__Kocuria"           "g__Pseudoalteromonas" "g__Capnocytophaga"    "g__Lonepinella"       "g__Butyrivibrio"
# otu_tab_common_freq <- otu_tab_common/(otu_tab_common + otu_tab_common[,"g__Prevotella"])
# shotgun_tab_common_freq <- shotgun_tab_common/(shotgun_tab_common + shotgun_tab_common[,"g__Prevotella"])


pdf(paste("plot_detected_genera_", var_name, "_rawRA_1_new.pdf", sep=""), height=7.7, width=7.7)
par(mfrow=c(4,4), pty="s", mar=c(2,3.5,2.5,0), cex.lab=1.3, font.lab=2)
for (i in 1:4) { # 1:4, 5:8
    # o <- order(var)
    # data.frame(otu_tab_common_freq[o,taxa_union[i]], shotgun_tab_common_freq[o,taxa_union[i]], var[o], otu_tab_common[o,taxa_union[i]], shotgun_tab_common[o,taxa_union[i]])
    
    main0 <- ifelse(i==1, "16S vs SMS RA", "")
    main1 <- ifelse(i==1, "16S RA", "")
    main2 <- ifelse(i==1, "SMS RA", "")
    main3 <- ifelse(i==1, "Average RA", "")
    xlab123 <- ifelse(i==6, var_name, "")
    
    plot(otu_tab_common_freq[,taxa_union[i]], shotgun_tab_common_freq[,taxa_union[i]], xlab="", ylab=unlist(strsplit(taxa_union[i], "_"))[3], main=main0)
    rho_value <- signif(cor(otu_tab_common_freq[,taxa_union[i]], shotgun_tab_common_freq[,taxa_union[i]]),2)
    legend(x="topright", legend=paste("cor=", rho_value))
    abline(coef(lm(shotgun_tab_common_freq[,taxa_union[i]] ~ otu_tab_common_freq[,taxa_union[i]])))
    abline(0,1,col="red")
    boxplot(otu_tab_common_freq[,taxa_union[i]] ~ var, main=main1, ylab="", xlab=xlab123)
    legend(x="topright", legend=paste("p =", ifelse(i %in% c(1,4,5,7), NA, signif(res.locom.16S$p.otu[,taxa_union[i]], 2))), bty="n")
    #legend(x="topright", legend=paste("pW =", signif(wilcox.test(otu_tab_common_freq[,taxa_union[i]] ~ var)$p.value, 2)), bty="n")
    boxplot(shotgun_tab_common_freq[,taxa_union[i]] ~ var, main=main2, ylab="", xlab=xlab123)
    legend(x="topright", legend=paste("p =", signif(res.locom.shotgun$p.otu[,taxa_union[i]], 2)), bty="n")
    #legend(x="topright", legend=paste("pW =", signif(wilcox.test(shotgun_tab_common_freq[,taxa_union[i]] ~ var)$p.value, 2)), bty="n")
    boxplot((otu_tab_common_freq[,taxa_union[i]] + shotgun_tab_common_freq[,taxa_union[i]])/2 ~ var, main=main3, ylab="", xlab=xlab123)
    legend(x="topright", legend=paste("p =", signif(res.locom.new$p.taxa.omni[,taxa_union[i]], 2)), bty="n")
    #legend(x="topright", legend=paste("pW =", signif(wilcox.test((otu_tab_common_freq[,taxa_union[i]] + shotgun_tab_common_freq[,taxa_union[i]])/2 ~ var)$p.value, 2)), bty="n")
}
dev.off()

# reference
ref <- match("g__Prevotella", colnames(otu_tab_common)) # 36
pdf(paste("plot_ref_genera_", var_name, ".pdf", sep=""), height=4, width=8)
par(mfrow=c(1,2), pty="s")

boxplot(otu_tab_common_freq[,ref] ~ var, xlab="", ylab="16S RA")
legend(x="topright", legend=paste("p =", signif(wilcox.test(otu_tab_common_freq[,ref] ~ var)$p.value, 2)), bty="n")
boxplot(shotgun_tab_common_freq[,ref] ~ var, xlab="", ylab="SMS RA")
legend(x="topright", legend=paste("p =", signif(wilcox.test(shotgun_tab_common_freq[,ref] ~ var)$p.value, 2)), bty="n")

dev.off()
