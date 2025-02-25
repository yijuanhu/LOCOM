#-----------------
# common taxa
#-----------------

common_taxa <- intersect(colnames(otu_tab_commonSam), colnames(shotgun_tab_commonSam))
length(common_taxa) # 236

# 16S 
otu_tab_common <- otu_tab_commonSam[, match(common_taxa, colnames(otu_tab_commonSam))]
dim(otu_tab_common) # 76 236
otu_tab_common[1:3, 1:3]
# g__Corynebacterium g__Coprococcus g__Streptococcus
# 11212.AC622                     0            262                6
# 11212.AC668                     3           1372              584
# 11212.Ahyb1000                  0            252                1
summary(rowSums(otu_tab_common))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 722    3064    4640    5281    7083   15491 

# shotgun
shotgun_tab_common <- shotgun_tab_commonSam[, match(common_taxa, colnames(shotgun_tab_commonSam))]
dim(shotgun_tab_common) # 76 236
shotgun_tab_common[1:3, 1:3]
# g__Corynebacterium g__Coprococcus g__Streptococcus
# 11212.AC622                     0            531              921
# 11212.AC668                     0            584            48644
# 11212.Ahyb1000                  8            389               17
summary(rowSums(shotgun_tab_common))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 6099   16828   33060  116604  140874  672637 

###################################
# var
###################################

var = as.factor(meta_common$diet)
var_name = "Diet"  
var_level = sort(unique(var))
var_nlevel = length(var_level)
table(var)
# var
# Folivore      Not
# 40       36

###################################
# Check overdispersion 
###################################

otu_tab_common_rrf <- Rarefy(otu_tab_common, depth = min(rowSums(otu_tab_common)))$otu.tab.rff # 722
shotgun_tab_common_rrf <- Rarefy(shotgun_tab_common, depth = min(rowSums(shotgun_tab_common)))$otu.tab.rff # 6032
discordance <- matrix(NA, nrow=length(common_sam), ncol=length(common_taxa))
discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf > 0] <- 0
discordance[otu_tab_common_rrf > 0 & shotgun_tab_common_rrf == 0] <- 1
discordance[otu_tab_common_rrf == 0 & shotgun_tab_common_rrf > 0] <- -1
discordance.col <- apply(discordance, 2, function(x) c(sum(x==1, na.rm=TRUE), sum(x==0, na.rm=TRUE), sum(x==-1, na.rm=TRUE))/sum(table(x)))
rowMeans(discordance.col, na.rm=TRUE)
#  0.1057703 0.1733124 0.7209174


###################################
# Scatter plots 
###################################

otu_tab_common_freq <- otu_tab_common/rowSums(otu_tab_common)
shotgun_tab_common_freq <- shotgun_tab_common/rowSums(shotgun_tab_common)
sort(colMeans(otu_tab_common_freq), decreasing = T)[1:5]
# g__Prevotella      g__Coprococcus g__Faecalibacterium     g__Ruminococcus      g__Clostridium 
# 0.17080156          0.10057741          0.07852142          0.07521036          0.06772252 
sort(colMeans(shotgun_tab_common_freq), decreasing = T)[1:5]
# g__Prevotella      g__Bacteroides      g__Clostridium     g__Ruminococcus g__Faecalibacterium 
# 0.16437445          0.11564247          0.05463376          0.05075046          0.05007498 
o <- order(colMeans(rbind(otu_tab_common_freq, shotgun_tab_common_freq)), decreasing = T)
genus_name_list <- colnames(otu_tab_common)

pdf("plot_scatter_Diet_1-15_36-50.pdf", height=18, width=15)
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
alpha_16S_r <- lm(alpha_16S~phylogency, data=meta_common)$residuals 
# shotgun
alpha_shotgun <- if (alpha_measure=="Shannon") diversity(shotgun_tab_common, index = "shannon", MARGIN = 1) else log(estimateR(shotgun_tab_common)[2,])
alpha_shotgun_r <- lm(alpha_shotgun~phylogency, data=meta_common)$residuals 
# comSam
alpha_comSam <- if (alpha_measure=="Shannon") diversity(tab_comSam, index = "shannon", MARGIN = 1) else log(estimateR(tab_comSam)[2,])
alpha_comSam_r <- lm(alpha_comSam~phylogency, data=meta_common_comSam)$residuals 
alpha_comSam_r <- matrix(alpha_comSam_r, ncol=1)
pvalue.comSam <- signif(ldm(alpha_comSam_r | strata ~ var_comSam, data=meta_common_comSam, seed=82955, dist.method="euclidean",
                            cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free",
                            scale.otu.table=FALSE, freq.scale.only=TRUE)$p.global.freq, 2)
# aveFreq
alpha_aveFreq <- if (alpha_measure=="Shannon") diversity(tab_aveFreq, index = "shannon", MARGIN = 1) else log(estimateR(tab_aveFreq)[2,])
alpha_aveFreq_r <- lm(alpha_aveFreq~phylogency, data=meta_common)$residuals 
# comCt
alpha_comCt <- if (alpha_measure=="Shannon") diversity(tab_comCt, index = "shannon", MARGIN = 1) else log(estimateR(tab_comCt)[2,])
alpha_comCt_r <- lm(alpha_comCt~phylogency, data=meta_common)$residuals 


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
# 0.96
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

p.perm.16S <- permanovaFL(otu_tab_common|phylogency ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.shotgun <- permanovaFL(shotgun_tab_common|phylogency ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.comSam <- permanovaFL(tab_comSam|phylogency + strata ~ var_comSam, data=meta_common_comSam, dist.method=dist_method, seed=82955,
                             cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free")$p.permanova
p.perm.aveFreq <- permanovaFL(tab_aveFreq|phylogency ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.comCt <- permanovaFL(tab_comCt|phylogency ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova

(p.com.p <- pcauchy(mean(c(tan( (0.5 - p.perm.16S)*pi), tan( (0.5 - p.perm.shotgun)*pi))), lower.tail=F))
# 0.00019996

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

dist_16S <- adjust.data.by.covariates(formula= ~ phylogency, data=meta_common, otu.table=otu_tab_common, dist.method=dist_method)$adj.dist 
dist_16S_eigen <- eigen(dist_16S, symmetric=TRUE)
dist_shotgun <- adjust.data.by.covariates(formula= ~ phylogency, data=meta_common, otu.table=shotgun_tab_common, dist.method=dist_method)$adj.dist 
dist_shotgun_eigen <- eigen(dist_shotgun, symmetric=TRUE)
dist_comSam <- adjust.data.by.covariates(formula= ~ phylogency, data=meta_common_comSam, otu.table=tab_comSam, dist.method=dist_method)$adj.dist 
dist_comSam_eigen <- eigen(dist_comSam, symmetric=TRUE)
dist_aveFreq <- adjust.data.by.covariates(formula= ~ phylogency, data=meta_common, otu.table=tab_aveFreq, dist.method=dist_method)$adj.dist 
dist_aveFreq_eigen <- eigen(dist_aveFreq, symmetric=TRUE)
dist_comCt <- adjust.data.by.covariates(formula= ~ phylogency, data=meta_common, otu.table=tab_comCt, dist.method=dist_method)$adj.dist 
dist_comCt_eigen <- eigen(dist_comCt, symmetric=TRUE)

dist_shotgun_eigen$vectors[,1] <- - dist_shotgun_eigen$vectors[,1]
dist_shotgun_eigen$vectors[,2] <- - dist_shotgun_eigen$vectors[,2]
dist_shotgun_eigen$vectors[,3] <- - dist_shotgun_eigen$vectors[,3]
dist_aveFreq_eigen$vectors[,1] <- - dist_aveFreq_eigen$vectors[,1]
dist_comCt_eigen$vectors[,1] <- - dist_comCt_eigen$vectors[,1]

pdf(paste("plot_beta_", dist_method, "_", var_name, ".pdf", sep=""), height=12, width=12)
par(mfrow=c(3,3), pty="s", mar=c(4,3,3,1))

xrange = c(-0.3, 0.34)
yrange = c(-0.3, 0.3)
gap = 0.015

# scatter plot of PC1
plot(dist_16S_eigen$vectors[,1], dist_shotgun_eigen$vectors[,1], main="PC1", xlab="16S data", ylab="SMS data", pch=16)
abline(0,1)
# scatter plot of PC2
plot(dist_16S_eigen$vectors[,2], dist_shotgun_eigen$vectors[,2], main="PC2", xlab="16S data", ylab="SMS data", pch=16)
abline(0,1)
# scatter plot of PC2
plot(dist_16S_eigen$vectors[,3], dist_shotgun_eigen$vectors[,3], main="PC3", xlab="16S data", ylab="SMS data", pch=16)
abline(0,1)

# 16S
plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="(a) 16S data", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
legend("topleft", legend=c("Folivore", "Not"), pch=c(1,16), col=c("blue","blue"))
legend("bottomright", legend=paste("p =", signif(p.perm.16S,3)), bty="n")

# shotgun
plot(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], main="(b) SMS data", xlab="PC1", ylab="PC2", col=color_shotgun, pch=pch_shotgun, xlim=xrange, ylim=yrange) # here use minus in y-axis to adjust 
text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
legend("topleft", legend=c("Folivore", "Not"), pch=c(1,16), col=c("red", "red"))
legend("bottomright", legend=paste("p =", signif(p.perm.shotgun,3)), bty="n")

# overlay
plot(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2], main="Overlaid figures (a) and (b)", xlab="PC1", ylab="PC2", col=color_16S, pch=pch_16S, xlim=xrange, ylim=yrange)
text(dist_16S_eigen$vectors[,1], dist_16S_eigen$vectors[,2]+gap, labels=ID, col=color_16S, cex=0.5)
ordiellipse(ord=dist_16S_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("blue", "blue"), lty=c(2,1))
points(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2], col=color_shotgun, pch=pch_shotgun) # here use minus in y-axis to adjust 
text(dist_shotgun_eigen$vectors[,1], dist_shotgun_eigen$vectors[,2]+gap, labels=ID, col=color_shotgun, cex=0.5)
ordiellipse(ord=dist_shotgun_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("red", "red"), lty=c(2,1))
legend("topleft", legend=c("16S, Folivore", "16S, Not", "SMS, Folivore", "SMS, Not"), 
       pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
legend("bottomright", legend=paste("p =", signif(p.com.p,3)), bty="n")

# combine by row
plot(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2], main="Combined observations", xlab="PC1", ylab="PC2", 
     col=c(color_16S, color_shotgun), pch=c(pch_16S, pch_shotgun), xlim=c(-0.2,0.2), ylim=c(-0.25,0.2))
text(dist_comSam_eigen$vectors[,1], dist_comSam_eigen$vectors[,2]+gap, labels=c(ID, ID), col=c(color_16S, color_shotgun), cex=0.5)
ordiellipse(ord=dist_comSam_eigen$vectors, groups=factor(var_comSam, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("bottomleft", legend=c("16S, Folivore", "16S, Not", "SMS, Folivore", "SMS, Not"), 
       pch=c(1,16,1,16), col=c("blue","blue","red","red"), cex=0.9)
legend("bottomright", legend=paste("p =", signif(p.perm.comSam,3)), bty="n")

# average frequencies
plot(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2], main="Averaged relative abundance", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
text(dist_aveFreq_eigen$vectors[,1], dist_aveFreq_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
ordiellipse(ord=dist_aveFreq_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("topleft", legend=c("Folivore", "Not"), pch=c(1,16), col=c("black", "black"))
legend("bottomright", legend=paste("p =", signif(p.perm.aveFreq,3)), bty="n")

# combine by pooling
plot(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2], main="Pooled count data", xlab="PC1", ylab="PC2", col=color_com, pch=pch_com, xlim=xrange, ylim=yrange)
text(dist_comCt_eigen$vectors[,1], dist_comCt_eigen$vectors[,2]+gap, labels=ID, col=color_com, cex=0.5)
ordiellipse(ord=dist_comCt_eigen$vectors, groups=factor(var, exclude=c()), conf=0.9, col=c("black", "black"), lty=c(2,1))
legend("topleft", legend=c("Folivore", "Not"), pch=c(1,16), col=c("black", "black"))
legend("bottomright", legend=paste("p =", signif(p.perm.comCt,3)), bty="n")

dev.off()


###################################
# LOCOM
###################################

source("../LOCOM_com_omni_fun.R")

otu_tab_common_filter <- otu_tab_common[, which(colMeans(otu_tab_common > 0) >= filter_thresh)] # 76 61
dim(otu_tab_common_filter) # 76 61
shotgun_tab_common_filter <- shotgun_tab_common[, which(colMeans(shotgun_tab_common > 0) >= filter_thresh)] 
dim(shotgun_tab_common_filter) # 76 189
var_num <- as.numeric(var)

# M2: dichotomous phylogency

meta_common$phylogency <- as.factor(meta_common$phylogency)
table(meta_common$phylogency)
# Apes   Lemurs NewWorld OldWorld 
# 9        9       29       29 
C <- as.factor(meta_common$phylogency)

fdr_target <- 0.1

#-------------
# 16S
#-------------
res.locom.16S <- locom(otu.table = otu_tab_common_filter, Y = var_num, C = C,
                       fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.16S$p.global # 9.999e-05
res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target] 
length(res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target]) 
# 24

length(res.locom.16S$q.otu) # 61

#-------------
# shotgun
#-------------
res.locom.shotgun <- locom(otu.table = shotgun_tab_common_filter, Y = var_num, C = C,
                           fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.shotgun$p.global # 0.00089991
res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target] 
length(res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target])
# 23
length(res.locom.shotgun$q.otu) # 189

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

length(q.comp[,which(q.comp < fdr_target)])  # 36
length(q.comp.HM[,which(q.comp.HM < fdr_target)])  # 36

# global p-value
pcauchy(mean(tan( (0.5-p.comp)*pi) ), lower.tail = F) # (10%) 0.002128551 
length(p.comp.HM)/sum(1/p.comp.HM) # (10%) 0.002107832 

length(p.comp) # 1062

#-------------
# Com-ct
#-------------
tab_comCt <- otu_tab_common + shotgun_tab_common
otus.keep.count <- which(colMeans(tab_comCt > 0) >= filter_thresh) 
length(otus.keep.count) # 195
tab_comCt_filter <- tab_comCt[, otus.keep.count]

res.locom.ct <- locom(otu.table = tab_comCt_filter, Y = var_num, C = C,
                      fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.ct$p.global # (10%) 0.00059994 
length(res.locom.ct$q.otu[,res.locom.ct$q.otu < fdr_target]) # 42 # 27 # 16 # 7

length(res.locom.ct$q.otu) # 195

#-------------
# Com-new
#-------------
length(otus.keep.count) # 195
mat1 <- 1:length(otus.keep.count)
mat2 <- 1:length(otus.keep.count)
res.locom.new <- Com2seq(taxa.table1 = otu_tab_common[,otus.keep.count], taxa.table2 = shotgun_tab_common[,otus.keep.count], Y1 = var_num, Y2 = var_num, C1 = C, C2 = C, mat1 = mat1, mat2 = mat2, 
                                cluster.id=factor(rep(1:n_sam, 2)), perm.within.type="none", perm.between.type="free",
                                fdr.nominal = fdr_target, seed=82955, n.perm.max=50000)
res.locom.new$p.global.wc # (10%) 9.999e-05 
length(res.locom.new$q.taxa.wc[,res.locom.new$q.taxa.wc < fdr_target])
# 30

res.locom.new$p.global.we # 9.999e-05
length(res.locom.new$q.taxa.we[,res.locom.new$q.taxa.we < fdr_target]) 
# 46

res.locom.new$p.global.omni # 9.999e-05 
length(res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < fdr_target]) 
# 54

length(res.locom.new$q.taxa.omni) # 195

save(res.locom.16S, res.locom.shotgun, res.locom.ct, res.locom.new, file="res.overlappingData.FDR10.RData")
load("res.overlappingData.FDR10.RData")


pdf(paste("plot_intercept_", var_name, ".pdf", sep=""), height=4, width=8)
par(mfrow=c(1,2), pty="s", mar=c(5,4,4,2))
plot(res.locom.new$para.est.wc[2,], res.locom.new$para.est.wc[3,], main="New-count", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16)
abline(0,1,col="red")
plot(res.locom.new$para.est.we[2,], res.locom.new$para.est.we[3,], main="New-equal", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16)
abline(0,1,col="red")
dev.off()


library(VennDiagram)
library(gridExtra)


set1 <- names(res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target])
# set1 <- names(res.locom.new$q.taxa.wc[,res.locom.new$q.taxa.wc < fdr_target])
set2 <- names(res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target])
# set2 <- names(res.locom.new$q.taxa.we[,res.locom.new$q.taxa.we < fdr_target])
# set3 <- names(res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < fdr_target])
set3 <- names(q.comp[,which(q.comp < fdr_target)])
# set3 <- names(res.locom.ct$q.otu[,res.locom.ct$q.otu < fdr_target])
label <- c("16S", "SMS", "Com-p-C")

grid.newpage()
venn.triple.plot <- draw.triple.venn(
    title = "FDR 10%",
    area1 = length(set1),
    area2 = length(set2),
    area3 = length(set3),
    n12 = length(intersect(set1, set2)),
    n13 = length(intersect(set1, set3)),
    n23 = length(intersect(set2, set3)),
    n123 = length(intersect(intersect(set1, set2), set3)),
    category = label,
    fill = "white",
    cat.col = "black",
    cat.cex = 3,
    cex = 3,
    cat.dist = c(0.1, 0.1, 0.1)
)

pdf("IBD_venn_3set_Com-p-C.pdf", width=9, height=9)
par(mfrow=c(1,1), pty="s")
grid.arrange(gTree(children=venn.triple.plot), top=textGrob("", gp=gpar(fontsize=60, font=1)))
dev.off()




