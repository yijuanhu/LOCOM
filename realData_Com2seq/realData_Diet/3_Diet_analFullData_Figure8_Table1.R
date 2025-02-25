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
# LOCOM
###################################

fdr_target <- 0.1

#-------------
# 16S
#-------------

Y <- ifelse(meta_otu$diet=="Folivore", 1, 0)
C <- as.factor(meta_otu$phylogency)
C <- model.matrix(~C)
res.locom.16S <- locom(otu.table = otu_tab_filter, Y = Y, C = C, 
                       fdr.nominal = fdr_target, seed=82955, n.perm.max=50000)
res.locom.16S$p.global # 0.000999001
res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target] 
length(res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target]) # 46 (after adjusting C), 70 (before adjusting C)
length(res.locom.16S$p.otu) # 122


res.locom.16S.C <- locom(otu.table = otu_tab_filter, Y = C, C = NULL, 
                       fdr.nominal = fdr_target, seed=82955, n.perm.max=50000)
res.locom.16S.C$p.global # 0.000999001
length(res.locom.16S.C$q.otu[,res.locom.16S.C$q.otu < fdr_target]) # 100

#-------------
# shotgun
#-------------

Y <- ifelse(meta_shotgun$diet=="Folivore", 1, 0)
C <- as.factor(meta_shotgun$phylogency)
C <- model.matrix(~C)
res.locom.shotgun <- locom(otu.table = shotgun_tab_filter, Y = Y, C = C, 
                           fdr.nominal = fdr_target, seed=82955, n.perm.max=50000) # 44000
res.locom.shotgun$p.global # 0.0003332223
res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target] 
length(res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target]) # 195 (after adjusting C), 246 (before adjusting C)
length(res.locom.shotgun$p.otu) # 1003

res.locom.shotgun.C <- locom(otu.table = shotgun_tab_filter, Y = C, C = NULL, 
                           fdr.nominal = fdr_target, seed=82955, n.perm.max=50000) # 44000
res.locom.shotgun.C$p.global # 0.000999001
length(res.locom.shotgun.C$q.otu[,res.locom.shotgun.C$q.otu < fdr_target]) # 476


#-------------
# Com-new
#-------------

Y1 <- ifelse(meta_otu$diet=="Folivore", 1, 0)
Y2 <- ifelse(meta_shotgun$diet=="Folivore", 1, 0)
C1 <- as.factor(meta_otu$phylogency)
C1 <- model.matrix(~C1)
C2 <- as.factor(meta_shotgun$phylogency)
C2 <- model.matrix(~C2)


# LOCOM package #
res.Com2seq <- Com2seq(table1 = otu_tab, table2 = shotgun_tab, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, 
                       seed = 82955, n.cores = 4, n.perm.max = 50000,
                       filter.thresh = 0.2, fdr.nominal = fdr_target)
res.Com2seq$p.global.wc # 9.999e-05
length(res.Com2seq$detected.taxa.wc) # 305

res.Com2seq$p.global.we # 9.999e-05
length(res.Com2seq$detected.taxa.we) # 335

res.Com2seq$p.global.omni # 0.00019998
length(res.Com2seq$detected.taxa.omni) # 382

dim(res.Com2seq$q.taxa.omni)
length(setdiff(res.Com2seq$detected.taxa.omni, res.Com2seq$detected.taxa.we)) # 53
length(setdiff(res.Com2seq$detected.taxa.omni, res.Com2seq$detected.taxa.wc)) # 80

length(res.locom.new$p.taxa.we) # 1062

# unadj
res.Com2seq.unadj <- Com2seq(table1 = otu_tab, table2 = shotgun_tab, Y1 = Y1, Y2 = Y2, C1 = NULL, C2 = NULL,
                       seed = 82955, n.cores = 4, n.perm.max = 50000,
                       filter.thresh = 0.2, fdr.nominal = fdr_target)
res.Com2seq.unadj$p.global.omni # 0.00019998
length(res.Com2seq.unadj$detected.taxa.omni) # 382


# confounder
res.Com2seq.C <- Com2seq(table1 = otu_tab, table2 = shotgun_tab, Y1 = C1, Y2 = C2, C1 = NULL, C2 = NULL, 
                       seed = 82955, n.cores = 4, n.perm.max = 50000,
                       filter.thresh = 0.2, fdr.nominal = fdr_target)
res.Com2seq.C$p.global.wc # 3.99992e-05
length(res.Com2seq.C$detected.taxa.wc) # 666

res.Com2seq.C$p.global.we # 3.99992e-05
length(res.Com2seq.C$detected.taxa.we) # 654

res.Com2seq.C$p.global.omni # 1.99996e-05
length(res.Com2seq.C$detected.taxa.omni) # 663 out of 1062



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

library(VennDiagram)
library(gridExtra)


set1 <- names(res.locom.16S$q.otu[,res.locom.16S$q.otu < fdr_target])
set2 <- names(res.locom.shotgun$q.otu[,res.locom.shotgun$q.otu < fdr_target])
set3 <- names(res.Com2seq$q.taxa.omni[,res.Com2seq$q.taxa.omni < fdr_target])
label <- c("16S", "SMS", "Com-2seq")

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

pdf("IBD_venn_3set_fullData_New-omni.pdf", width=9, height=9)
par(mfrow=c(1,1), pty="s")
grid.arrange(gTree(children=venn.triple.plot), top=textGrob("", gp=gpar(fontsize=60, font=1)))
dev.off()







