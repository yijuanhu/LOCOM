% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LOCOM_fun.R
\name{locom}
\alias{locom}
\title{A logistic regression model for testing differential abundance in compositional microbiome data (LOCOM)}
\usage{
locom(
  otu.table,
  Y,
  C = NULL,
  fdr.nominal = 0.2,
  seed = NULL,
  n.perm.max = NULL,
  n.rej.stop = 100,
  n.cores = 4,
  filter.thresh = 0.2,
  permute = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{otu.table}{the OTU table (or taxa count table) in which the rows correspond to samples and the columns correspond to OTUs (taxa).}

\item{Y}{the trait of interest. It can be an array, matrix, or data frame; a factor should be represented by the corresponding design matrix. All components in a matrix or data frame will be tested for microbial association jointly.}

\item{C}{the other (confounding) covariates to be adjusted. See other requirements in \code{Y}.}

\item{fdr.nominal}{the nominal FDR. The default is 0.2.}

\item{seed}{a user-supplied integer seed for the random number generator in the
permutation procedure. The default is NULL, which means that an integer seed will be
generated internally and randomly. In either case, the integer seed will be stored
in the output object in case the user wants to reproduce the permutation replicates.}

\item{n.perm.max}{the maximum number of permutations. The default is NULL, in which case \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal})
is used where \code{n.otu} is the total number of OTUs (that have non-zero counts in at least one sample).}

\item{n.rej.stop}{the minimum number of rejections (i.e., the permutation test
statistic exceeds the observed test statistic) to obtain before stopping the permutation procedure. The
default is 100.}

\item{n.cores}{the number of cores to be used for parallel computing. The default is 1.}

\item{filter.thresh}{a real value between 0 and 1; OTUs present in fewer than \code{filter.thresh} * 100\% samples are filtered out. The default is 0.2.}

\item{permute}{a logical value indicating whether to perform permutation. The default is TRUE.}

\item{verbose}{a logical value indicating whether to generate verbose output during the permutation process. Default is TRUE.}
}
\value{
A list consisting of
\itemize{
  \item effect.size - effect size at each OTU, i.e., beta_j,1 - median_j'=1,...J(beta_j',1)
  \item p.otu - p-values for OTU-specific tests
  \item q.otu - q-values (adjusted p-values by the BH procedure) for OTU-specific tests
  \item detected.otu - detected OTUs (using the column names of the OTU table) at the nominal FDR based on \code{q.otu}
  \item p.global - p-value for the global test
  \item n.perm.completed - number of permutations completed
  \item seed - the seed used to generate the permutation replicates
  \item p.otu.asymptotic - asymptotic p-values for OTU-specific tests
  \item q.otu.asymptotic - asymptotic q-values (adjusted p-values by the BH procedure) for OTU-specific tests
  \item detected.otu.asymptotic - detected OTUs (using the column names of the OTU table) at the nominal FDR based on \code{q.otu.asymptotic}
}
}
\description{
This function allows you to test
(1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on the log ratio of relative abundances between pairs of taxa, and
(2). whehter the whole community is associated with the trait (a global test).
The tests accommodate continuous, discrete (binary, categorical), and multivariate traits, and allows adjustment of confounders.
}
\details{
This function uses a sequential stopping criterion (similar to that of Sandve et al. 2011) for the permutation procedure,
which stops when all taxon-level tests have either reached the pre-specified
number (default 100) of rejections or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) that is below the
nominal FDR level (default 0.2). The permutation procedure is always terminated if a pre-specified maximum number (see description of \code{n.perm.max}) of
permutations have been generated. The global test is based on all permutation replicates when the procedure stops/terminates.
}
\examples{
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
filter.out.otu <- which(prop.presence < 0.4) # In practice, use threshold 0.2; 
                                             # here, filter out more OTUs for speed
if (length(filter.out.otu) > 0) {
    throat.otu.table.filter <- throat.otu.table.filter[, -filter.out.otu]
    prop.presence <- prop.presence[-filter.out.otu]
}

# running locom
res <- locom(otu.table = throat.otu.table.filter, Y = Y, C = C, 
             seed = 123, n.cores = 1, n.perm.max = 1000)
res$p.global 
res$detected.otu
}
