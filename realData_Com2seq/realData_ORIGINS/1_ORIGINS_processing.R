#########################################################################################################
# ORIGINS: investigate the association between periodontal bacteria and prediabetes prevalence
#########################################################################################################

library(vtable)
library(vegan)
library(GUniFrac)
library(readr)
library(biomformat)
library(prodlim)
library(parallel)
library(matrixStats)
library(LOCOM)

filter_thresh <- 0.2


#------------------------
# meta data
#------------------------

meta <- read_delim("ID11808/11808_20220510-095139.txt")
meta <- data.frame(meta)
dim(meta) # 384  67
table(meta$prediabetes) # Prediabetes is the trait in this project
# 0              1   Not applicable   not provided 
#259             57             62              6 
sample_sub <- which(meta$prediabetes==0 | meta$prediabetes==1)
meta <- meta[sample_sub,]
o1 <- order(meta$sample_name)
meta <- meta[o1,]
dim(meta) # 316 67
table(meta$prediabetes)
# 0   1 
# 259  57
meta$prediabetes[meta$prediabetes=="1"] = "Prediabetes"
meta$prediabetes[meta$prediabetes=="0"] = "Not"
table(meta$prediabetes)
# Not Prediabetes 
# 259          57 
# non-informative/redundant variables
# var_tab <- st(meta, out = "return") # check variable information
meta <- meta[, -c(8, 13, 15, 16, 18:26, 31:44, 46, 47, 53, 55, 56, 58, 59, 61, 62, 64, 65)]
dim(meta) # 316  29
colnames(meta)
# "aapperio"    "aapperio4" (4 categories):?   
# "ahei1" "ahei2" "ahei3" "ahei4": Maybe the AHEI grades your diet, assigning a score ranging from 0 (nonadherence) to 110 (perfect adherence), based on how often you eat certain foods, both healthy and unhealthy fare.     
# "apdqs" "apdqst"(categorical based on apdqs(0-55; 55-66; >65)): Maybe is the APDQS score, ranging from 0 to 132, was calculated from increasing intake of 20 beneficial components (including fruit, vegetables, legumes, low-fat dairy, fish, and moderate alcohol intake) and decreasing intake of 13 adverse components (including fried foods, salty snacks, desserts, high-fat dairy, and sugar-sweetened soft drinks).     
# "bmi" "height" "weightkg"        
# "cigcurr": smoking status     
# "crp": C-creative protein
# "educ": education        
# "glucosecrc": plasma glucose  
# "hba1c": hemoglobin A1C      
# "host_age": age    
# "insulinv1": maybe insulin
# "meanaloss": mean attachment loss   
# "meandbp" "meansbp": mean blood pressure systolic and diastolic 
# "meanpd": mean probing depth      
# "mets": metabolic equivalent       
# "perbop": ?     
# "prediabetes"
# "raceethn": race/ethnicity  
# "sex"         
# "totmisst": tooth loss    

#------------------------
# 16S data
#------------------------

otu_biom <- read_biom("ID11808/141016_otu_table.biom")
otu_tab <- t(as.matrix(biom_data(otu_biom)))
dim(otu_tab) # 366 1853

# genus names
otu_taxa <- observation_metadata(otu_biom)
# group to a higher taxonomic level 
level <- 6 # 6--genus, 5--family, 4--order, 3--class, 2--phylum, 1--kindom
unique_taxa <- unique(otu_taxa[,1:level])  
dim(unique_taxa) # 348   6
mat <- row.match(otu_taxa[,1:level], unique_taxa)
otu_tab1 <- NULL
for (i in 1:nrow(unique_taxa)) {
  otu_tab1 <- cbind(otu_tab1, rowSums(otu_tab[,which(mat==i), drop=FALSE]))
}
dim(otu_tab1) # 366 348

otu_tab <- otu_tab1
tax_tab <- unique_taxa
tax_tab[which(tax_tab[,6]=="g__"),6] <- tax_tab[which(tax_tab[,6]=="g__"),5]
tax_tab[which(tax_tab[,5]=="f__"),6] <- tax_tab[which(tax_tab[,5]=="f__"),4]
tax_tab[which(tax_tab[,4]=="o__"),6] <- tax_tab[which(tax_tab[,4]=="o__"),3]
tax_tab[which(tax_tab[,3]=="c__"),6] <- tax_tab[which(tax_tab[,3]=="c__"),2]
colnames(otu_tab) <- tax_tab[,6]
dim(otu_tab) # 366 348
#otu_tab[,"g__Ignavigranum"]
otu_tab[,"g__Chelonobacter"]
otu_tab[1:3,1:3]
#f__Sphingobacteriaceae g__Porphyromonas g__Akkermansia
#11808.4587.SALB.W2V1                      1              212              4
#11808.3441.SALB.W2V1                      0              678              0
#11808.3542.SALB.W2V1                      0             1130              0

# remove samples with low library sizes
summary(rowSums(otu_tab))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1    8407   23517   20393   29691   49842
sort(rowSums(otu_tab))[1:100]
w <- which(rowSums(otu_tab) < 5000)
otu_tab <- otu_tab[-w, colSums(otu_tab)>=5]
dim(otu_tab) # 276 234

# match with meta data
int1 <- intersect(meta$sample_name, rownames(otu_tab))
length(int1) # 271
mat11 <- match(int1, meta$sample_name)
meta_otu <- meta[mat11,]
mat12 <- match(int1, rownames(otu_tab))
otu_tab <- otu_tab[mat12, ]
dim(meta_otu) # 271  29
dim(otu_tab) # 271 234
summary(rowSums(otu_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8041   21728   27188   26950   31232   49842 


# filter rare taxa
otu_tab_filter <- otu_tab[, which(colMeans(otu_tab > 0) >= filter_thresh)]
dim(otu_tab_filter) # 271  85
summary(rowSums(otu_tab_filter))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8041   21710   27181   26922   31210   49809 
sort(colMeans(otu_tab_filter/rowSums(otu_tab_filter)) , decreasing = T)[1:5]
# g__Haemophilus    g__Prevotella g__Streptococcus   g__Veillonella   f__Neisseriaceae 
# 0.26236882       0.20376384       0.10612820       0.09356362       0.07956955 

# # permanova
# res.perm.16S.filter <- permanovaFL(otu_tab_filter ~ prediabetes, data=meta_otu, dist.method="bray", seed=82955)
# res.perm.16S.filter$p.permanova

#------------------------
# shotgun data
#------------------------

shotgun_biom <- read_biom("ID11808/120040_genus.biom")
shotgun_tab <- t(as.matrix(biom_data(shotgun_biom)))
dim(shotgun_tab) # 192 1091

# genus names
taxa_name_split <- strsplit(colnames(shotgun_tab),";")
genus_shotgun <- sapply(taxa_name_split, tail, 1)
unique_genus_shotgun <- unique(genus_shotgun)
length(unique_genus_shotgun) # 1091, all unique!
colnames(shotgun_tab) <- genus_shotgun
shotgun_tab[1:3,1:3]
# g__Methanobrevibacter g__Acidiplasma g__Haloplasma
# 11808.3160.SALB.W2V1                     0              0             0
# 11808.3070.SALB.W2V1                     0              0             0
# 11808.3069.SALB.W2V1                     0              0             0

# remove samples with low library sizes
summary(rowSums(shotgun_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 25   45744  132538  168126  250484 1161306 
w <- which(rowSums(shotgun_tab) < 5000)
shotgun_tab  <- shotgun_tab[-w, colSums(shotgun_tab)>=5]
dim(shotgun_tab) # 183 756

# match with meta data
int2 <- intersect(meta$sample_name, rownames(shotgun_tab))
length(int2) # 183
mat21 <- match(int2, meta$sample_name)
meta_shotgun <- meta[mat21,]
mat22 <- match(int2, rownames(shotgun_tab))
shotgun_tab <- shotgun_tab[mat22, ]
dim(meta_shotgun) # 183  29
dim(shotgun_tab) # 183 756
summary(rowSums(shotgun_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8890   56950  143192  176321  255426 1161240 

# filter rare taxa
shotgun_tab_filter <- shotgun_tab[, which(colMeans(shotgun_tab > 0) >= filter_thresh)]
dim(shotgun_tab_filter) # 183 403
summary(rowSums(shotgun_tab_filter))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8890   56942  143169  176258  255176 1160842 
sort(colMeans(shotgun_tab_filter/rowSums(shotgun_tab_filter)) , decreasing = T)[1:5]
# g__Prevotella g__Streptococcus     g__Neisseria   g__Haemophilus  g__Veillonella 
# 0.23784059       0.13799035       0.12832197       0.12281151       0.08698881

# res.perm.shotgun.filter <- permanovaFL(shotgun_tab_filter ~ prediabetes, data=meta_shotgun, dist.method="bray", seed=82955)
# res.perm.shotgun.filter$p.permanova


##########################################################
# Common samples, common taxa
##########################################################

#-----------------
# common sample
#-----------------

common_sam <- intersect(rownames(otu_tab), rownames(shotgun_tab))
length(common_sam) # 152

# meta data
meta_common <- meta[match(common_sam, meta$sample_name),]
table(meta_common$prediabetes)
# 0   1 
# 105  47 

# 16S data
otu_tab_commonSam <- otu_tab[match(common_sam, rownames(otu_tab)), ] 
dim(otu_tab_commonSam) # 152 234

# shotgun data
shotgun_tab_commonSam <- shotgun_tab[match(common_sam, rownames(shotgun_tab)), ] 
dim(shotgun_tab_commonSam) # 152 756

# samples unique to 16S and shotgun
meta_otu_unique <- meta_otu[-match(common_sam, meta_otu$sample_name),] 
dim(meta_otu_unique) # 119 29
table(meta_otu_unique$prediabetes)
# Not Prediabetes 
# 118           1 

meta_shotgun_unique <- meta_shotgun[-match(common_sam, meta_shotgun$sample_name),] 
dim(meta_shotgun_unique) # 31 29
table(meta_shotgun_unique$prediabetes)
# Not Prediabetes 
# 23           8 

152+119+31
# 302

length(union(rownames(otu_tab), rownames(shotgun_tab)))
# 302
length(union(colnames(otu_tab), colnames(shotgun_tab)))
# 864


# meta_union

meta_union <- rbind(meta_common, meta_otu_unique, meta_shotgun_unique)
dim(meta_union) # 302 29
table(meta_union$prediabetes)
# Not Prediabetes 
# 246          56 
w0 <- (meta_union$prediabetes=="Not")
w1 <- (meta_union$prediabetes=="Prediabetes")

# age, sex, race/ethnicity, education, smoking status, body mass index, and baseline glucose levels
age <- as.numeric(meta_union$host_age)
wilcox.test(age[w0], age[w1]) # 8.86e-06***

chisq.test(meta_union$sex, meta_union$prediabetes) # 0.1476

chisq.test(meta_union$raceethn, meta_union$prediabetes) # 0.0171***

chisq.test(meta_union$educ, meta_union$prediabetes) # 0.2299

table(meta_union$cigcurr)
chisq.test(meta_union$cigcurr, meta_union$prediabetes) # 0.1652

bmi <- as.numeric(meta_union$bmi)
wilcox.test(bmi[w0], bmi[w1]) # 8.86e-06***

glu <- as.numeric(meta_union$glucosecrc)
wilcox.test(glu[w0], glu[w1]) # 4.163e-08***

insulin <- as.numeric(meta_union$insulinv1)
wilcox.test(insulin[w0], insulin[w1]) # 0.002597



