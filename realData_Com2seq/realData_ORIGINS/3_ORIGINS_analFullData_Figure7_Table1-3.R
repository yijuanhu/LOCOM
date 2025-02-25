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
# combining data
###################################

n_sam <- nrow(meta_common)
ID <- 1:n_sam
ID_comSam <- c(ID, ID)
var_comSam <- as.factor(c(var, var))
meta_common_comSam <- rbind(meta_common, meta_common)
meta_common_comSam$strata <- c(rep(0, n_sam), rep(1, n_sam))

# tab_comSam <- rbind(otu_tab_common, shotgun_tab_common)

###################################
# Alpha diversity
###################################

alpha_measure <- "Shannon"  # Shannon, log Chao1

# 16S
alpha_16S <- if (alpha_measure=="Shannon") diversity(otu_tab_commonSam, index = "shannon", MARGIN = 1) else log(estimateR(otu_tab_commonSam)[2,])
alpha_16S_r <- lm(alpha_16S~0, data=meta_common)$residuals 
# shotgun
alpha_shotgun <- if (alpha_measure=="Shannon") diversity(shotgun_tab_commonSam, index = "shannon", MARGIN = 1) else log(estimateR(shotgun_tab_commonSam)[2,])
alpha_shotgun_r <- lm(alpha_shotgun~0, data=meta_common)$residuals 
# comSam
alpha_comSam <- c(alpha_16S, alpha_shotgun)
alpha_comSam_r <- lm(alpha_comSam~0, data=meta_common_comSam)$residuals 
alpha_comSam_r <- matrix(alpha_comSam_r, ncol=1)
pvalue.comSam <- signif(ldm(alpha_comSam_r | strata ~ var_comSam, data=meta_common_comSam, seed=82955, dist.method="euclidean",
                            cluster.id = ID_comSam, perm.within.type = "none", perm.between.type = "free",
                            scale.otu.table=FALSE, freq.scale.only=TRUE)$p.global.freq, 2)

pdf(paste("plot_alpha_", alpha_measure, "_", var_name, "_fullTaxa.pdf", sep=""), height=6.5, width=6.5)
par(mfrow=c(2,2), pty="s", mar=c(4.5,3,4,0.1))

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
# 0.086
# comSam
boxplot(alpha_comSam_r ~ var_comSam, main="Combined observations", xlab="", ylab=alpha_measure) 
legend(x="bottomright", legend=paste("p =", pvalue.comSam), bty="n")

dev.off()


###################################
# Beta diversity
###################################t

dist_method = "bray"

p.perm.16S <- permanovaFL(otu_tab_commonSam ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
p.perm.shotgun <- permanovaFL(shotgun_tab_commonSam ~ var, data=meta_common, dist.method=dist_method, seed=82955)$p.permanova
(p.com.p <- pcauchy(mean(c(tan( (0.5 - p.perm.16S)*pi), tan( (0.5 - p.perm.shotgun)*pi))), lower.tail=F))
# 0.117105

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

dist_16S <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=otu_tab_commonSam, dist.method=dist_method)$adj.dist 
dist_16S_eigen <- eigen(dist_16S, symmetric=TRUE)
dist_shotgun <- adjust.data.by.covariates(formula= ~ 1, data=meta_common, otu.table=shotgun_tab_commonSam, dist.method=dist_method)$adj.dist 
dist_shotgun_eigen <- eigen(dist_shotgun, symmetric=TRUE)

dist_16S_eigen$vectors[,2] <- - dist_16S_eigen$vectors[,2]
dist_shotgun_eigen$vectors[,2] <- - dist_shotgun_eigen$vectors[,2]
dist_shotgun_eigen$vectors[,3] <- - dist_shotgun_eigen$vectors[,3]

pdf(paste("plot_beta_", dist_method, "_", var_name, "_fullTaxa.pdf", sep=""), height=12, width=12)
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

dev.off()


###################################
# LOCOM
###################################

source("../LOCOM_com_omni_fun_general.R")

# M1: no confounder
C <- NULL

fdr_target <- 0.1

#-------------
# 16S
#-------------

Y <- ifelse(meta_otu$prediabetes=="Prediabetes", 1, 0)
res.locom.16S <- locom(otu.table = otu_tab_filter, Y = Y, C = C, 
                       fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000)
res.locom.16S$p.global # (20%) 0.06824431 (10%) 0.06824431
res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target] 
# g__Chelonobacter 
# 0.06374469 

# g__Chelonobacter 
# 0.06374469 

length(res.locom.16S$p.otu) # 85

#-------------
# shotgun
#-------------

Y <- ifelse(meta_shotgun$prediabetes=="Prediabetes", 1, 0)
res.locom.shotgun <- locom(otu.table = shotgun_tab_filter, Y = Y, C = C, 
                           fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.perm.max=50000) # 50000
res.locom.shotgun$p.global # 0.08804345  0.08661827
res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target] 
# g__Gemella g__Ignavigranum 
# 0.1328038       0.1328038 

# numeric(0)

length(res.locom.shotgun$p.otu) # 403

#-------------
# Com-p 
#-------------
name1 <- colnames(res.locom.16S$p.otu)
name2 <- colnames(res.locom.shotgun$p.otu)
common.name <- intersect(name1, name2)
j.mat1 <- match(common.name, name1)
j.mat2 <- match(common.name, name2)

p.both <- rbind(res.locom.16S$p.otu[j.mat1], res.locom.shotgun$p.otu[j.mat2])
p.comp <- matrix(pcauchy(apply(tan((0.5-p.both)*pi), 2, mean), lower.tail = F), nrow = 1)
p.comp.HM <- 2/colSums(1/p.both)
p.comp.name <- common.name
if (length(res.locom.16S$p.otu[-j.mat1])>0) {
    p.comp <- c(p.comp, res.locom.16S$p.otu[-j.mat1])
    p.comp.HM <- c(p.comp.HM, res.locom.16S$p.otu[-j.mat1])
    p.comp.name <-c(p.comp.name, name1[-j.mat1])
}
if (length(res.locom.shotgun$p.otu[-j.mat2])>0) {
    p.comp <- c(p.comp, res.locom.shotgun$p.otu[-j.mat2])
    p.comp.HM <- c(p.comp.HM, res.locom.shotgun$p.otu[-j.mat2])
    p.comp.name <-c(p.comp.name, name2[-j.mat2])
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

q.comp[,which(q.comp < fdr_target)]  
# g__Chelonobacter g__Pseudarthrobacter           g__Gemella      g__Ignavigranum 
# 0.1874957            0.1874957            0.1449967            0.1449967 

# numeric(0)

q.comp.HM[,which(q.comp.HM < fdr_target)]  
# g__Chelonobacter g__Pseudarthrobacter           g__Gemella      g__Ignavigranum 
# 0.1874957            0.1874957            0.1449967            0.1449967 

# numeric(0)

# global p-value
pcauchy(mean(tan( (0.5-p.comp)*pi) ), lower.tail = F) # 0.2620908  0.2456133
length(p.comp.HM)/sum(1/p.comp.HM) # 0.0467922  0.04667769

length(q.comp) # 440

#-------------
# Com-ct
#-------------

sam.union <- union(rownames(otu_tab), rownames(shotgun_tab))
n.sam.union <- length(sam.union)
taxa.union <- union(colnames(otu_tab), colnames(shotgun_tab))
n.taxa.union <- length(taxa.union)

otu_tab_fill0 <- matrix(0, nrow = n.sam.union, ncol = n.taxa.union)
shotgun_tab_fill0 <- matrix(0, nrow = n.sam.union, ncol = n.taxa.union)

rownames(otu_tab_fill0) <- sam.union
rownames(shotgun_tab_fill0) <- sam.union
colnames(otu_tab_fill0) <- taxa.union
colnames(shotgun_tab_fill0) <- taxa.union

otu_tab_fill0[match(rownames(otu_tab), sam.union), match(colnames(otu_tab), taxa.union)] <- otu_tab
shotgun_tab_fill0[match(rownames(shotgun_tab), sam.union), match(colnames(shotgun_tab), taxa.union)] <- shotgun_tab
pooled_tab <- otu_tab_fill0 + shotgun_tab_fill0

pooled_tab_filter <- pooled_tab[, which(colMeans(pooled_tab > 0) >= filter_thresh)]
dim(pooled_tab_filter) # 302 365

Y <- rep(NA, n.sam.union)
Y[match(rownames(otu_tab), sam.union)] <- meta_otu$prediabetes
Y[match(rownames(shotgun_tab), sam.union)] <- meta_shotgun$prediabetes

Y <- ifelse(Y=="Prediabetes", 1, 0)
res.locom.ct <- locom(otu.table = pooled_tab_filter, Y = Y, C = C,
                      fdr.nominal = fdr_target, weight.type = "count", seed=82955, n.cores = 1, n.perm.max=50000)

res.locom.ct$p.global # 0.00039996  0.00039996
res.locom.ct$detected.otu
# [1] "g__Akkermansia"           "f__Lachnospiraceae"       "f__Neisseriaceae"         "g__Catonella"            
# [5] "g__[Prevotella]"          "g__Leptotrichia"          "o__Bacteroidales"         "g__Moryella"             
# [9] "o__Clostridiales"         "p__SR1"                   "f__Peptostreptococcaceae" "g__Peptococcus"          
# [13] "f__Planococcaceae"        "f__Cardiobacteriaceae"    "f__[Weeksellaceae]"       "g__N09"                  
# [17] "f__Actinomycetaceae"      "g__TG5"                   "f__Gemellaceae"           "f__Comamonadaceae"       
# [21] "f__[Mogibacteriaceae]"    "g__Schwartzia"            "g__Geobacillus"           "g__Campylobacter"        
# [25] "o__Burkholderiales"       "f__Pasteurellaceae"       "f__Carnobacteriaceae"     "g__Butyrivibrio"         
# [29] "g__Thermus"               "g__Bulleidia"             "g__Paludibacter"          "o__RF39"                 
# [33] "f__Aerococcaceae"         "g__Gemella"               "g__Arthrobacter"          "g__Anaerovorax"          
# [37] "g__Listeria"              "g__Serinibacter"          "g__Curtobacterium"        "g__Acaricomes"           
# [41] "g__Nesterenkonia"         "g__Pseudarthrobacter"     "g__Zhihengliuella"        "g__Novibacillus"         
# [45] "g__Ignavigranum"          "g__Weissella"             "g__Clostridiisalibacter"  "g__Negativicoccus"       
# [49] "g__Salmonella"  

# [1] "g__Akkermansia"           "f__Lachnospiraceae"       "f__Neisseriaceae"         "g__Catonella"            
# [5] "g__[Prevotella]"          "g__Leptotrichia"          "o__Bacteroidales"         "g__Moryella"             
# [9] "o__Clostridiales"         "p__SR1"                   "f__Peptostreptococcaceae" "g__Peptococcus"          
# [13] "f__Planococcaceae"        "f__[Weeksellaceae]"       "g__N09"                   "f__Actinomycetaceae"     
# [17] "f__Gemellaceae"           "f__Comamonadaceae"        "f__[Mogibacteriaceae]"    "g__Schwartzia"           
# [21] "g__Geobacillus"           "g__Campylobacter"         "o__Burkholderiales"       "f__Pasteurellaceae"      
# [25] "f__Carnobacteriaceae"     "g__Butyrivibrio"          "g__Paludibacter"          "f__Aerococcaceae"        
# [29] "g__Gemella"               "g__Anaerovorax"           "g__Listeria"              "g__Serinibacter"         
# [33] "g__Curtobacterium"        "g__Acaricomes"            "g__Pseudarthrobacter"     "g__Zhihengliuella"       
# [37] "g__Novibacillus"          "g__Ignavigranum"          "g__Clostridiisalibacter"  "g__Negativicoccus"

length(res.locom.ct$detected.otu) # 49  40

#-------------
# Com-new
#-------------

Y1 <- ifelse(meta_otu$prediabetes=="Prediabetes", 1, 0)
Y2 <- ifelse(meta_shotgun$prediabetes=="Prediabetes", 1, 0)
res.locom.new <-Com2seq(table1 = otu_tab, table2 = shotgun_tab, Y1 = Y1, Y2 = Y2, C1 = C, C2 = C, 
                               filter_thresh = 0.2,
                               fdr.nominal = fdr_target,
                               n.perm.max = 50000, n.rej.stop = 100, n.cores = 1, ref.taxon = NULL, seed=82955)
length(res.locom.new$p.taxa.wc) # 440
res.locom.new$p.global.wc # 0.02129787  0.03875922
res.locom.new$q.taxa.wc[,res.locom.new$q.taxa.wc < fdr_target] 
# [1] "g__Campylobacter"        "g__Butyrivibrio"         "g__Gemella"              "f__Lachnospiraceae"     
# [5] "g__Moryella"             "f__Comamonadaceae"       "g__Actinotignum"         "g__Serinibacter"        
# [9] "g__Curtobacterium"       "g__Acaricomes"           "g__Nesterenkonia"        "g__Pseudarthrobacter"   
# [13] "g__Psychromicrobium"     "g__Renibacterium"        "g__Zhihengliuella"       "g__Novibacillus"        
# [17] "g__Ignavigranum"         "g__Weissella"            "g__Clostridiisalibacter" "g__Negativicoccus"      
# [21] "g__Salmonella"   

# g__Gemella g__Ignavigranum 
# 0.05719886      0.05719886 
length(res.locom.new$detected.taxa.wc) # 21  2

res.locom.new$p.global.we # 0.07619238   0.02623948
res.locom.new$q.taxa.we[,res.locom.new$q.taxa.we < fdr_target] 
# [1] "g__Butyrivibrio"         "f__Lachnospiraceae"      "g__Moryella"             "f__Comamonadaceae"      
# [5] "g__Actinotignum"         "g__Serinibacter"         "g__Curtobacterium"       "g__Acaricomes"          
# [9] "g__Nesterenkonia"        "g__Pseudarthrobacter"    "g__Psychromicrobium"     "g__Renibacterium"       
# [13] "g__Zhihengliuella"       "g__Novibacillus"         "g__Ignavigranum"         "g__Weissella"           
# [17] "g__Clostridiisalibacter" "g__Negativicoccus"       "g__Salmonella"   

# g__Butyrivibrio g__Ignavigranum 
# 0.04839903      0.04399912 
length(res.locom.new$detected.taxa.we) # 19  2

res.locom.new$p.global.omni # 0.04389561 0.03397932
res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < fdr_target] 
# g__Campylobacter         g__Butyrivibrio              g__Gemella      f__Lachnospiraceae             g__Moryella 
# 0.21371001              0.10695171              0.09679806              0.10695171              0.10695171 
# f__Comamonadaceae         g__Actinotignum         g__Serinibacter       g__Curtobacterium           g__Acaricomes 
# 0.17831222              0.17831222              0.10695171              0.10695171              0.10695171 
# g__Nesterenkonia    g__Pseudarthrobacter     g__Psychromicrobium        g__Renibacterium       g__Zhihengliuella 
# 0.19975600              0.10695171              0.10695171              0.17831222              0.13434398 
# g__Novibacillus         g__Ignavigranum            g__Weissella g__Clostridiisalibacter       g__Negativicoccus 
# 0.10695171              0.07919842              0.17831222              0.10695171              0.10695171 
# g__Salmonella  
# 0.12948312 

# g__Butyrivibrio      g__Gemella g__Ignavigranum 
# 0.07626514      0.07626514      0.06159877 
length(res.locom.new$detected.taxa.omni) # 21  3

res.locom.new$q.taxa.omni[,res.locom.new$q.taxa.omni < 0.3] 

res.locom.new$q.taxa.omni[,c("g__Pseudoalteromonas", "g__Actinomyces", "g__Capnocytophaga", "g__Lonepinella", "g__Campylobacter", "g__Kocuria")]

length(res.locom.new$q.taxa.omni) # 440

# all(res.locom.new$detected.taxa.omni==res.locom.new$detected.taxa.wc)
# TRUE
# setdiff(res.locom.new$detected.taxa.omni, res.locom.new$detected.taxa.we)
# "g__Campylobacter" "g__Gemella" 

save(res.locom.16S, res.locom.shotgun, res.locom.ct, res.locom.new, file="res.fullData.FDR10.RData")
load("res.fullData.FDR10.RData")









pdf(paste("plot_intercept_", var_name, "_fulldata.pdf", sep=""), height=4, width=8)
par(mfrow=c(1,2), pty="s", mar=c(5,4,4,2))
plot(res.locom.new$para.est.wc[2,], res.locom.new$para.est.wc[3,], main="New-wc", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16)
abline(0,1,col="red")
plot(res.locom.new$para.est.we[2,], res.locom.new$para.est.we[3,], main="New-we", 
     xlab="Intercept for 16S data", ylab="Intercept for SMS data", pch=16)
abline(0,1,col="red")
dev.off()


#################
# summarizing
#################

taxa_union <- c("g__Butyrivibrio", "g__Gemella", "g__Ignavigranum", "g__Chelonobacter")
taxa_union_16S <- c("g__Butyrivibrio", "g__Chelonobacter")
signif(res.locom.new$p.taxa.omni[, taxa_union],3)
signif(res.locom.new$p.taxa.we[, taxa_union],3)
signif(res.locom.new$p.taxa.wc[, taxa_union],3)
signif(p.comp[, taxa_union],3)
signif(p.comp.HM[, taxa_union],3)
signif(res.locom.16S$p.otu[, taxa_union_16S],3)
signif(res.locom.shotgun$p.otu[, taxa_union],3)

signif(res.locom.new$q.taxa.omni[, taxa_union],3)
signif(res.locom.new$q.taxa.we[, taxa_union],3)
signif(res.locom.new$q.taxa.wc[, taxa_union],3)
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

otu_tab <- data.frame(otu_tab, g__Ignavigranum=0)
otu_tab_freq <- otu_tab/rowSums(otu_tab)
shotgun_tab_freq <-shotgun_tab/rowSums(shotgun_tab)
    
pdf(paste("plot_detected_genera_", var_name, "_rawRA_fullData.pdf", sep=""), height=7.7, width=3.9)
par(mfrow=c(4,2), pty="s", mar=c(2,3.5,2.5,0), cex.lab=1.3, font.lab=2)
for (i in 1:4) { 
    main1 <- ifelse(i==1, "16S RA", "")
    main2 <- ifelse(i==1, "SMS RA", "")
    
    boxplot(otu_tab_freq[,taxa_union[i]] ~ meta_otu$prediabetes, main=main1, ylab=unlist(strsplit(taxa_union[i], "_"))[3], xlab=xlab123)
    legend(x="topright", legend=paste("p =", ifelse(i %in% c(2,3), NA, signif(res.locom.16S$p.otu[,taxa_union[i]], 2))), bty="n")
    #legend(x="topright", legend=paste("pW =", signif(wilcox.test(otu_tab_common_freq[,taxa_union[i]] ~ var)$p.value, 2)), bty="n")
    boxplot(shotgun_tab_freq[,taxa_union[i]] ~ meta_shotgun$prediabetes, main=main2, ylab="", xlab=xlab123)
    legend(x="topright", legend=paste("p =", signif(res.locom.shotgun$p.otu[,taxa_union[i]], 2)), bty="n")
    #legend(x="topright", legend=paste("pW =", signif(wilcox.test(shotgun_tab_freq[,taxa_union[i]] ~ meta_shotgun$prediabetes)$p.value, 2)), bty="n")
}
dev.off()

