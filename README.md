# LOCOM

A logistic regression model for testing differential abundance in compositional microbiome data

The LOCOM package provides tests of differential abundance at individual taxa (or OTUs) and in a microbial community as a whole, that is based on the log ratio of relative abundances between each taxon and a reference taxon. While the 'locom' function was developed to analyze one sequencing dataset, the newer 'Com2seq' function combines two sequencing datasets (e.g., from 16S marker-gene and shotgun metagenomic sequencing) for integrative analysis of differential abundance. All tests accommodate continuous, discrete (binary, categorical), and multivariate traits, and allows adjustment of confounders.

## Installation

```r
devtools::install_github("yijuanhu/LOCOM")
```

## An example of using the locom function

Apply 'locom' to a dataset from the study of smoking effects on the microbiome of the human upper respiratory tract (Charlson et al., 2010):

```r
library(LOCOM)
data("throat.otu.table")
data("throat.meta")
data("throat.otu.taxonomy")

# filtering out three samples with antibiotic use
filter.out.sam <- which(throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")
throat.otu.table.filter <- throat.otu.table[-filter.out.sam,]
throat.meta.filter <- throat.meta[-filter.out.sam,]
Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)

# filtering out rare OTUs
prop.presence <- colMeans(throat.otu.table.filter > 0)
filter.out.otu <- which(prop.presence < 0.2)
if (length(filter.out.otu) > 0) {
    throat.otu.table.filter <- throat.otu.table.filter[, -filter.out.otu]
    prop.presence <- prop.presence[-filter.out.otu]
}

# running locom
res <- locom(otu.table = throat.otu.table.filter, Y = Y, C = C, seed = 1, n.cores=4, n.perm.max = 1000)
res$p.global
res$detected.otu

# summarizing results
w <- match(res$detected.otu, colnames(res$p.otu))
o <- w[order(res$p.otu[w])]

summary.table <- data.frame(otu.name = colnames(res$p.otu)[o],
                            mean.freq = colMeans(throat.otu.table.filter/rowSums(throat.otu.table.filter))[o],
                            prop.presence = prop.presence[o],
                            p.value = signif(res$p.otu[o], 3),
                            q.value = signif(res$q.otu[o], 3),
                            effect.size = signif(res$effect.size[o], 3),
                            otu.tax = throat.otu.taxonomy[as.numeric(colnames(res$p.otu)[o]) + 1],
                            row.names = NULL)
summary.table
```

## An example of using the Com2seq function

Apply 'Com2seq' to a simulated dataset:

```r
library(LOCOM)
data("otu.table.1")
data("otu.table.2")
data("meta.1")
data("meta.2")

# preparing Y1, Y2, C1, C2
Y1 <- cbind(1*(meta.1[,"Y"]==1), 1*(meta.1[,"Y"]==2)) # get design matrix for a three-level factor
Y2 <- cbind(1*(meta.2[,"Y"]==1), 1*(meta.2[,"Y"]==2))
C1 <- meta.1[,"C"] 
C2 <- meta.2[,"C"] 

# restricting to a subset of samples and OTUs for speed
sub.sam.1 <- c(31:60)
sub.sam.2 <- c(1:20, 51:60)
otu.table.1 <- otu.table.1[sub.sam.1, 1:70]
otu.table.2 <- otu.table.2[sub.sam.2, 1:70]
Y1 <- Y1[sub.sam.1,]
C1 <- C1[sub.sam.1]
Y2 <- Y2[sub.sam.2,]
C2 <- C2[sub.sam.2]

# running Com2seq
res.Com2seq <- Com2seq(table1 = otu.table.1, table2 = otu.table.2, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, 
                       seed = 123, n.cores = 1, n.perm.max = 1000)
res.Com2seq$p.global.omni 
res.Com2seq$detected.taxa.omni 
```
