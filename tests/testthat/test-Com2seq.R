context("Testing `Com2seq` function")
library(LOCOM)
library(testthat)

data("otu.table.1")
data("otu.table.2")
data("meta.1")
data("meta.2")

# preparing Y1, Y2, C1, C2
Y1 <- model.matrix(~meta.1[,"Y"])[,-1] # get design matrix for a three-level factor
Y2 <- model.matrix(~meta.2[,"Y"])[,-1]
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

# test
test_that("`Com2seq` function provides expected results", {
    res.Com2seq <- Com2seq(table1 = otu.table.1, table2 = otu.table.2, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2, 
                           seed = 123, n.cores = 1, n.perm.max = 1000)
    res_p <- signif(res.Com2seq$p.global.omni, 3)
    expect_equivalent(res_p, 0.049)
})
