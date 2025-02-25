#########################################################################################################
# Evolutionary trends in host physiology outweigh dietary niche in structuring primate gut microbiomes
#########################################################################################################
 
library(vtable)
library(vegan)
library(GUniFrac)
library(readr)
library(biomformat)
library(prodlim)
library(parallel)
library(matrixStats)
library(gtsummary)

library(LOCOM)

filter_thresh <- 0.2

#------------------------
# meta data
#------------------------

meta <- read_delim("ID11212/11212_20170918-134804.txt")
meta <- data.frame(meta)
dim(meta) # 172  41
# var_tab <- st(meta, out = "return")
meta <- meta[, -c(2:4, 5:9, 11, 13, 15:18, 21, 22, 25, 28, 30:41)]
dim(meta) # 172  11
colnames(meta)
# "sample_name" "latitude"
# "sample_name"          "country"           "diet"
# "elevation"            "geo_loc_name"
# "habitat2"             "host_common_name"
# "host_scientific_name" "host_taxid"
# "latitude"             "longitude"

# tbl <- tbl_summary(meta[,-1]) 
# as_gt(tbl) %>% gt::as_latex()
# as_kable(tbl, format = "latex")
table(meta$diet) # Diet is the trait in this project
# Folivore      Not 
# 94       78
o <- order(meta$sample_name)
meta <- meta[o,]
dim(meta) # 172  11

Apes_w <- which(meta$host_scientific_name=="Pan troglodytes"
                |meta$host_scientific_name=="Gorilla gorilla")
length(Apes_w)

Lemurs_w <- which(meta$host_scientific_name=="Eulemur rubriventer"
                  |meta$host_scientific_name=="Lemur catta"
                  |meta$host_scientific_name=="Propithecus verreauxi")
length(Lemurs_w)

OldWorld_w <- which(meta$host_scientific_name=="Theropithecus gelada"
                    |meta$host_scientific_name=="Papio anubis"
                    |meta$host_scientific_name=="Papio hamadryas"
                    |meta$host_scientific_name=="Cercopithecus ascanius"
                    |meta$host_scientific_name=="Colobus guereza"
                    |meta$host_scientific_name=="Piliocolobus badius")
length(OldWorld_w)

meta$phylogency <- rep("NewWorld", nrow(meta))
meta$phylogency[Apes_w] <-"Apes"
meta$phylogency[Lemurs_w] <-"Lemurs"
meta$phylogency[OldWorld_w] <- "OldWorld"
table(meta$phylogency)
# Apes   Lemurs NewWorld OldWorld 
# 24       33       70       45 
meta$phylogency.bin <- meta$phylogency
meta$phylogency.bin[which(meta$phylogency.bin!="NewWorld")] <- "Else"

#------------------------
# 16S data
#------------------------

otu_biom <- read_biom("ID11212/53182_otu_table.biom")
otu_tab <- t(as.matrix(biom_data(otu_biom)))
dim(otu_tab) # 154 9464

# genus names
otu_taxa <- observation_metadata(otu_biom)
# group to a higher taxonomic level 
level <- 6 # 6--genus, 5--family, 4--order, 3--class, 2--phylum, 1--kindom
unique_taxa <- unique(otu_taxa[,1:level])  
dim(unique_taxa) # 841  6
mat <- row.match(otu_taxa[,1:level], unique_taxa)
otu_tab1 <- NULL
for (i in 1:nrow(unique_taxa)) {
  otu_tab1 <- cbind(otu_tab1, rowSums(otu_tab[,which(mat==i), drop=FALSE]))
}
dim(otu_tab1) # 154 841

otu_tab <- otu_tab1
rm(otu_biom, otu_tab1)
tax_tab <- unique_taxa
tax_tab[which(tax_tab[,6]=="g__"),6] <- tax_tab[which(tax_tab[,6]=="g__"),5]
tax_tab[which(tax_tab[,5]=="f__"),6] <- tax_tab[which(tax_tab[,5]=="f__"),4]
tax_tab[which(tax_tab[,4]=="o__"),6] <- tax_tab[which(tax_tab[,4]=="o__"),3]
tax_tab[which(tax_tab[,3]=="c__"),6] <- tax_tab[which(tax_tab[,3]=="c__"),2]
colnames(otu_tab) <- tax_tab[,6]
dim(otu_tab) # 154 841
otu_tab[1:3,1:3]
# o__Bacteroidales g__Corynebacterium g__Sharpea
# 11212.Gor1              1168                  0          0
# 11212.Gor11             2210                  0          0
# 11212.Gor12             2457                  0          0             

# remove samples with low library sizes
summary(rowSums(otu_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4262   11528   16646   16992   21865   36021 
sort(rowSums(otu_tab))[1:100]
w <- which(rowSums(otu_tab) < 5000)
otu_tab <- otu_tab[-w, colSums(otu_tab)>=5]
dim(otu_tab) # 153 524

# match with meta data
int1 <- intersect(meta$sample_name, rownames(otu_tab))
length(int1) # 153
mat11 <- match(int1, meta$sample_name)
meta_otu <- meta[mat11,]
mat12 <- match(int1, rownames(otu_tab))
otu_tab <- otu_tab[mat12, ]
dim(meta_otu) # 153  13
dim(otu_tab) # 153 524
summary(rowSums(otu_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6121   11602   16670   17072   21878   36011 

# filter rare taxa
otu_tab_filter <- otu_tab[, which(colMeans(otu_tab > 0) >= filter_thresh)]
dim(otu_tab_filter) # 153 122
summary(rowSums(otu_tab_filter))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6079   11489   16437   16765   21347   35796 
sort(colMeans(otu_tab_filter/rowSums(otu_tab_filter)) , decreasing = T)[1:5]
# o__Clostridiales f__Ruminococcaceae f__Lachnospiraceae            o__RF39      g__Prevotella 
# 0.18680473         0.16669453         0.07083562         0.05027651         0.04946008 

# # permanova
# res.perm.16S.filter <- permanovaFL(otu_tab_filter ~ diet, data=meta_otu, dist.method="bray", seed=82955)
# res.perm.16S.filter$p.permanova # 0.00019996

#------------------------
# shotgun data
#------------------------

shotgun_biom <- read_biom("ID11212/115839_genus.biom")
shotgun_tab <- t(as.matrix(biom_data(shotgun_biom)))
dim(shotgun_tab) # 95 2254

# genus names
taxa_name_split <- strsplit(colnames(shotgun_tab),";")
genus_shotgun <- sapply(taxa_name_split, tail, 1)
unique_genus_shotgun <- unique(genus_shotgun)
length(unique_genus_shotgun) # 2254, all unique!
colnames(shotgun_tab) <- genus_shotgun
shotgun_tab[1:3,1:3]
# g__Natrialba g__Natronococcus g__Natronolimnobius
# 11212.APAL36            0                0                   0
# 11212.AC668             0                0                   0
# 11212.APAL38            0                0                   0
summary(rowSums(shotgun_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 9957   35195   58064  167111  194840 2409072
rm(shotgun_biom, taxa_name_split)
shotgun_tab  <- shotgun_tab[-w, colSums(shotgun_tab)>=5]
dim(shotgun_tab) # 95 1778

# match with meta data
int2 <- intersect(meta$sample_name, rownames(shotgun_tab))
length(int2) # 95
mat21 <- match(int2, meta$sample_name)
meta_shotgun <- meta[mat21,]
mat22 <- match(int2, rownames(shotgun_tab))
shotgun_tab <- shotgun_tab[mat22, ]
dim(meta_shotgun) # 95 13
dim(shotgun_tab) # 95 1778
summary(rowSums(shotgun_tab))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9951   35191   58058  167100  194836 2409026 


# filter rare taxa
shotgun_tab_filter <- shotgun_tab[, which(colMeans(shotgun_tab > 0) >= filter_thresh)]
dim(shotgun_tab_filter) # 95 1003
summary(rowSums(shotgun_tab_filter))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 9901   35114   57979  166476  194696 2369932
sort(colMeans(shotgun_tab_filter/rowSums(shotgun_tab_filter)) , decreasing = T)[1:5]
# g__Prevotella      g__Bacteroides      g__Clostridium g__Faecalibacterium        g__Treponema 
# 0.09938050          0.07441563          0.04264408          0.03277757          0.03065230 

# res.perm.shotgun.filter <- permanovaFL(shotgun_tab_filter ~ diet, data=meta_shotgun, dist.method="bray", seed=82955)
# res.perm.shotgun.filter$p.permanova # 0.00019996

##########################################################
# Common samples, common taxa
##########################################################

#-----------------
# common sample
#-----------------

common_sam <- intersect(rownames(otu_tab), rownames(shotgun_tab))
length(common_sam) # 76

# meta data
meta_common <- meta[match(common_sam, meta$sample_name),]
table(meta_common$diet)
# Folivore      Not 
# 40       36

# 16S data
otu_tab_commonSam <- otu_tab[match(common_sam, rownames(otu_tab)), ] 
dim(otu_tab_commonSam) # 76 524

# shotgun data
shotgun_tab_commonSam <- shotgun_tab[match(common_sam, rownames(shotgun_tab)), ] 
dim(shotgun_tab_commonSam) # 76 1778

# samples unique to 16S and shotgun
meta_otu_unique <- meta_otu[-match(common_sam, meta_otu$sample_name),] 
dim(meta_otu_unique) # 77 13
table(meta_otu_unique$diet)
# Folivore      Not 
# 45       32 

meta_shotgun_unique <- meta_shotgun[-match(common_sam, meta_shotgun$sample_name),] 
dim(meta_shotgun_unique) # 19 13
table(meta_shotgun_unique$diet)
# Folivore      Not 
# 9       10 

length(union(rownames(otu_tab), rownames(shotgun_tab)))
# 172
length(union(colnames(otu_tab), colnames(shotgun_tab)))
# 2062


# meta_union

meta_union <- rbind(meta_common[,c("diet", "phylogency")], meta_otu_unique[,c("diet", "phylogency")], meta_shotgun_unique[,c("diet", "phylogency")])
dim(meta_union) # 172 2
chisq.test(meta_union$diet, meta_union$phylogency) # p-value=1.091e-07

