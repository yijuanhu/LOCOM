context("Testing `locom` function")
library(LOCOM)
library(testthat)

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
filter.out.otu <- which(prop.presence < 0.4) # In practice, use threshold 0.2; here, filter out more OTUs for speed
if (length(filter.out.otu) > 0) {
    throat.otu.table.filter <- throat.otu.table.filter[, -filter.out.otu]
    prop.presence <- prop.presence[-filter.out.otu]
}

# test
test_that("`locom` function provides expected results", {
    res <- locom(otu.table = throat.otu.table.filter, Y = Y, C = C, seed = 123, n.cores = 1, n.perm.max = 1000)
    res_p <- signif(res$p.global, 3)
    expect_equivalent(res_p, 0.011)
})
