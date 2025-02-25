#' Com2seq: combining two sequencing datasets (e.g., from 16S marker-gene and shotgun metagenomic sequencing) for integrative analysis of differential abundance at both the genus and community levels
#'
#' This function allows you to test
#' (1). whether any genus is associated with the trait of interest with FDR control, based on the log ratio of relative abundances between pairs of taxa, and
#' (2). whether the whole community is associated with the trait (a global test).
#' The tests accommodate continuous, discrete (binary, categorical), and multivariate traits, and allows adjustment of confounders.
#'
#' This function uses a sequential stopping criterion (similar to that of Sandve et al. 2011) for the permutation procedure,
#' which stops when all taxon-level tests have either reached the pre-specified
#' number (default 100) of rejections or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) that is below the
#' nominal FDR level (default 0.2). The permutation procedure is always terminated if a pre-specified maximum number (see description of \code{n.perm.max}) of
#' permutations have been generated. The global test is based on all permutation replicates when the procedure stops/terminates.
#'
#'
#' @param table1 taxa count (relative abundance) table agglomerated at the genus level. The rows correspond to samples and the columns correspond to genera.
#' @param table2 taxa count (relative abundance) table agglomerated at the genus level. The rows correspond to samples and the columns correspond to genera.
#' @param Y1 trait of interest, with samples ordered in the same way as in \code{table1}. It can be an array, matrix, or data frame; a factor should be represented by the corresponding design matrix. All components in a matrix or data frame will be tested for microbial association jointly.
#' @param Y2 trait of interest, with samples ordered in the same way as in \code{table2}. See other requirements in \code{Y1}.
#' @param C1 other covariates, with samples ordered in the same way as in \code{table1}. See other requirements in \code{Y1}.
#' @param C2 other covariates, with samples ordered in the same way as in \code{table2}. See other requirements in \code{Y1}.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param seed a user-supplied integer seed for the random number generator in the
#'   permutation procedure. The default is NULL, which means that an integer seed will be
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case the user wants to reproduce the permutation replicates.
#' @param n.perm.max the maximum number of permutations. The default is NULL, in which case \code{n.taxa} * \code{n.rej.stop} * (1/\code{fdr.nominal})
#'   is used where \code{n.taxa} is the total number of taxa (that have non-zero counts in at least one sample).
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation test
#'   statistic exceeds the observed test statistic) to obtain before stopping the permutation procedure. The
#'   default is 100.
#' @param n.cores the number of cores to be used for parallel computing. The default is 1.
#' @param iter_max Maximum number of iterations for the Newton algorithm. The default is 200.
#' @param clean A logical value indicating whether to remove taxa that fail to converge. Default is FALSE.
#'   When set to TRUE, taxa that do not converge will be excluded from the analysis.
#' @param filter.thresh a real value between 0 and 1; OTUs present in fewer than \code{filter.thresh} * 100\% samples are filtered out. The default is 0.2.
#' @param verbose a logical value indicating whether to generate verbose output during the permutation process. Default is TRUE.
#' @return A list consisting of
#' \itemize{
#'   \item effect.size.wc - effect size at each taxon, i.e., beta_j,1 - median_j'=1,...J (beta_j',1), based on the count weight
#'   \item p.taxa.wc - p-values by the taxon-specific tests based on the count weight
#'   \item q.taxa.wc - q-values (adjusted p-values by the BH procedure) by the taxon-specific tests based on the count weight
#'   \item detected.taxa.wc - detected taxa (using the column names of the taxon table) at the nominal FDR by the taxon-specific tests based on the count weight
#'   \item p.global.wc - p-value by the global test based on the count weight
#'   \item effect.size.we - effect size at each taxon, i.e., beta_j,1 - median_j'=1,...J(beta_j',1), based on the equal weight
#'   \item p.taxa.we - p-values by the taxon-specific tests based on the equal weight
#'   \item q.taxa.we - q-values (adjusted p-values by the BH procedure) by the taxon-specific tests based on the equal weight
#'   \item detected.taxa.we - detected taxa (using the column names of the taxon table) at the nominal FDR by the taxon-specific tests based on the equal weight
#'   \item p.global.we - p-value by the global test based on the equal weight
#'   \item p.taxa.omni - p-values by the taxon-specific omnibus tests
#'   \item q.taxa.omni - q-values (adjusted p-values by the BH procedure) by the taxon-specific omnibus tests
#'   \item detected.taxa.omni - detected taxa (using the column names of the taxon table) at the nominal FDR by the taxon-specific omnibus tests
#'   \item p.global.omni - p-value by the global omnibus test
#'   \item n.perm.completed - number of permutations completed
#'   \item ref.taxon - the selected reference taxon
#'   \item seed - the seed used to generate the permutation replicates
#' }
#' @export
#' @examples
#' library(LOCOM)
#' data("otu.table.1")
#' data("otu.table.2")
#' data("meta.1")
#' data("meta.2")
#'
#' # preparing Y1, Y2, C1, C2
#' Y1 <- cbind(1*(meta.1[,"Y"]==1), 1*(meta.1[,"Y"]==2)) # get design matrix for a three-level factor
#' Y2 <- cbind(1*(meta.2[,"Y"]==1), 1*(meta.2[,"Y"]==2))
#' C1 <- meta.1[,"C"]
#' C2 <- meta.2[,"C"]
#'
#' # restricting to a subset of samples and OTUs for speed
#' sub.sam.1 <- c(31:60)
#' sub.sam.2 <- c(1:20, 51:60)
#' otu.table.1 <- otu.table.1[sub.sam.1, 1:70]
#' otu.table.2 <- otu.table.2[sub.sam.2, 1:70]
#' Y1 <- Y1[sub.sam.1,]
#' C1 <- C1[sub.sam.1]
#' Y2 <- Y2[sub.sam.2,]
#' C2 <- C2[sub.sam.2]
#'
#' # running Com2seq
#' res.Com2seq <- Com2seq(table1 = otu.table.1, table2 = otu.table.2,
#'                        Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2,
#'                        seed = 123, n.cores = 1, n.perm.max = 1000)
#' res.Com2seq$p.global.omni
#' res.Com2seq$detected.taxa.omni

Com2seq.cor <- function(table1, table2, Y1, Y2, C1 = NULL, C2 = NULL,
                    fdr.nominal = 0.2, seed = NULL,
                    n.perm.max = NULL, n.rej.stop = 100, n.cores = 1,
                    iter_max = 200, clean = FALSE,
                    filter.thresh = 0.2, verbose = TRUE) {
  
  Firth.thresh = 0.4
  tol = 1e-6
  robust_var = FALSE
  we = FALSE
  
  if (is.data.frame(table1)) table1 <- data.matrix(table1)
  if (is.data.frame(table2)) table2 <- data.matrix(table2)
  if (is.data.frame(Y1)) Y1 <- data.matrix(Y1)
  if (is.data.frame(Y2)) Y2 <- data.matrix(Y2)
  if (is.data.frame(C1)) C1 <- data.matrix(C1)
  if (is.data.frame(C2)) C2 <- data.matrix(C2)
  
  n.trait <- length(Y1)/nrow(table1)
  Y1 <- matrix(Y1, ncol=n.trait)
  Y2 <- matrix(Y2, ncol=n.trait)
  if (!is.null(C1)) {
    n.conf <- length(C1)/nrow(table1)
    C1 <- matrix(C1, ncol=n.conf)
    C2 <- matrix(C2, ncol=n.conf)
  }
  
  # detect whether the data are relative abundances
  # if the data are relative abundances, then we just apply the equal weights across all data
  if (mean(rowSums(table1)) < 10 | mean(rowSums(table2)) < 10) {
    we.only = TRUE
  } else {
    we.only = FALSE
  }
  
  # sort samples in the two tables
  
  n.sam1 <- nrow(table1)
  n.sam2 <- nrow(table2)
  sam.common <- intersect(rownames(table1), rownames(table2))
  n.sam.common <- length(sam.common)
  
  n.sam1.only <- n.sam1 - n.sam.common
  n.sam2.only <- n.sam2 - n.sam.common
  
  sam1.sort <- sam.common
  sam2.sort <- sam.common
  if (n.sam1.only > 0) sam1.sort <- c(setdiff(rownames(table1), sam.common), sam.common)
  if (n.sam2.only > 0) sam2.sort <- c(sam.common, setdiff(rownames(table2), sam.common))
  
  mat1 <- match(sam1.sort, rownames(table1)) # indices
  mat2 <- match(sam2.sort, rownames(table2)) # indices
  
  table1 <- table1[mat1, ]
  table2 <- table2[mat2, ]
  
  Y1 <- Y1[mat1,,drop=FALSE]
  Y2 <- Y2[mat2,,drop=FALSE]
  Y <- rbind(Y1, Y2)
  
  if (is.null(C1) | is.null(C2)) {
    C <- NULL
  } else {
    C1 <- C1[mat1,,drop=FALSE]
    C2 <- C2[mat2,,drop=FALSE]
    C <- rbind(C1, C2)
  }
  
  # create three tables (table.comTaxa.filter, table1.only.filter, table2.only.filter)
  
  taxa.common <- intersect(colnames(table1), colnames(table2))
  
  table1.filter <- table1[ ,which(colMeans(table1 > 0) >= filter.thresh)] # filter 1
  table2.filter <- table2[ ,which(colMeans(table2 > 0) >= filter.thresh)] # filter 2
  
  taxa1.only.filter <- setdiff(colnames(table1.filter), taxa.common) # only taxa in taxa1
  taxa2.only.filter <- setdiff(colnames(table2.filter), taxa.common) # only taxa in taxa2
  n.taxa1.only.filter <- length(taxa1.only.filter)
  n.taxa2.only.filter <- length(taxa2.only.filter)
  if (n.taxa1.only.filter > 0) {
    table1.only.filter <- table1.filter[, taxa1.only.filter, drop=FALSE]
  }
  if (n.taxa2.only.filter > 0) {
    table2.only.filter <- table2.filter[, taxa2.only.filter, drop=FALSE]
  }
  rm(taxa1.only.filter, taxa2.only.filter)
  
  sam.union <- union(rownames(table1), rownames(table2))
  table1.uniSam.comTaxa <- matrix(0, nrow = length(sam.union), ncol = length(taxa.common))
  table2.uniSam.comTaxa <- matrix(0, nrow = length(sam.union), ncol = length(taxa.common))
  table1.uniSam.comTaxa[match(rownames(table1), sam.union), ] <- table1[, match(taxa.common, colnames(table1))] # fill the table1's rows elements
  table2.uniSam.comTaxa[match(rownames(table2), sam.union), ] <- table2[, match(taxa.common, colnames(table2))]
  table.uniSam.comTaxa <- table1.uniSam.comTaxa + table2.uniSam.comTaxa
  
  taxa.common.filter3 <- taxa.common[which(colMeans(table.uniSam.comTaxa > 0) >= filter.thresh)] # filter 3
  
  rm(table1.uniSam.comTaxa, table2.uniSam.comTaxa, table.uniSam.comTaxa, sam.union)
  
  taxa.3filters <- unique(c(colnames(table1.filter), colnames(table2.filter), taxa.common.filter3)) # passing three filters
  taxa.common.filter <- intersect(taxa.common, taxa.3filters)
  
  table.comTaxa.filter <- rbind(table1[, match(taxa.common.filter, colnames(table1))],
                                table2[, match(taxa.common.filter, colnames(table2))]) # same samples and same taxon
  
  n.sam.com <- nrow(table.comTaxa.filter) # samples of two source
  n.taxa.common.filter <- ncol(table.comTaxa.filter)
  
  # find reference OTU
  ref.taxon <- NULL
  if (is.null(ref.taxon)) {
    freq.table.comTaxa.filter <- table.comTaxa.filter/rowSums(table.comTaxa.filter)
    mean.freq <- colMeans(freq.table.comTaxa.filter)
    sd.freq <- colSds(freq.table.comTaxa.filter)
    ref.taxon <- which.max(mean.freq - sd.freq)
    rm(freq.table.comTaxa.filter)
  }
  table.comTaxa.filter[table.comTaxa.filter[, ref.taxon]==0, ref.taxon] <- 1
  
  # -----------------------
  # Observed statistic
  # -----------------------
  
  int1 <- c(rep(1, n.sam1), rep(0, n.sam2))
  int2 <- c(rep(0, n.sam1), rep(1, n.sam2))
  
  if (is.null(C) == TRUE) {
    X <- cbind(Y, int1, int2) # int1 and int2 is the intercept (in the supp, they are the first 2 elements)
    Yr <- resid(lm(Y ~ -1 + int1 + int2))
    Yr <- matrix(Yr, ncol=length(Yr)/nrow(Y))
    
    X1 <- cbind(Y1, 1)
    Y1r <- Y1
    
    X2 <- cbind(Y2, 1)
    Y2r <- Y2
  } else {
    X <- cbind(Y, int1, int2, C)
    Yr <- resid(lm(Y ~ -1 + int1 + int2 + C))
    Yr <- matrix(Yr, ncol=length(Yr)/nrow(Y))
    
    X1 <- cbind(Y1, 1, C1)
    Y1r <- resid(lm(Y1 ~ C1))
    Y1r <- matrix(Y1r, ncol=length(Y1r)/nrow(Y1))
    
    X2 <- cbind(Y2, 1, C2)
    Y2r <- resid(lm(Y2 ~ C2))
    Y2r <- matrix(Y2r, ncol=length(Y2r)/nrow(Y2))
  }
  
  freq.table.ref <- table.comTaxa.filter/(table.comTaxa.filter + table.comTaxa.filter[, ref.taxon]) # Y_j/(Y_j+Y_J)
  N_denominator = table.comTaxa.filter + table.comTaxa.filter[, ref.taxon] # Y_j+Y_J #weight.count
  dim(freq.table.ref)
  dim(N_denominator)
  
  n.X <- ncol(X)
  zero.vector = rep(0,(ncol(X)-2))
  row_vector1 = c(zero.vector,c(-1,1))
  row_vector2 = c(zero.vector,c(1,-1))
  # I1 = matrix(rep(row_vector1, n.sam1), nrow = nrow(X1), byrow = TRUE)
  # I2 = matrix(rep(row_vector2, n.sam2), nrow = nrow(X2), byrow = TRUE)
  
  if (!we.only) weight.count = (table.comTaxa.filter + table.comTaxa.filter[, ref.taxon])
  if(we) weight.equal = matrix(rep(colMeans(table.comTaxa.filter + table.comTaxa.filter[, ref.taxon]), n.sam.com), ncol=n.taxa.common.filter, byrow=TRUE) # temporary
  prop.presence <- colMeans(table.comTaxa.filter > 0) # rarity analysis
  
  if (n.taxa1.only.filter > 0) {
    freq.table1.ref <- table1.only.filter/(table1.only.filter + table.comTaxa.filter[1:n.sam1, ref.taxon])
    N_denominator1 = table1.only.filter + table.comTaxa.filter[1:n.sam1, ref.taxon]
    n.X1 <- ncol(X1)
    #XX1 <- CalculateXY(X1,X1)
    weight1 = (table1.only.filter + table.comTaxa.filter[1:n.sam1, ref.taxon])
    prop.presence1 <- colMeans(table1.only.filter>0)
  }
  if (n.taxa2.only.filter > 0) {
    freq.table2.ref <- table2.only.filter/(table2.only.filter + table.comTaxa.filter[-(1:n.sam1), ref.taxon])
    N_denominator2 = table2.only.filter + table.comTaxa.filter[-(1:n.sam1), ref.taxon]
    n.X2 <- ncol(X2)
    #XX2 <- CalculateXY(X2,X1)
    weight2 = (table2.only.filter + table.comTaxa.filter[-(1:n.sam1), ref.taxon])
    prop.presence2 <- colMeans(table2.only.filter>0)
  }
  beta_init0 <- array(0, dim = c(n.X, n.taxa.common.filter))
  
  # get the estimated rho
  dim(freq.table.ref)
  dim(N_denominator)
  dim(beta_init0)
  dim(weight.count)
  beta_init0_all = beta_init0
  beta_init0_all[,ref.taxon] = c(0.00001,0.00001,0.00001)
  res.obs.wc.all <- Newton_cor(freq_table = freq.table.ref, N_denominator = N_denominator, X = X,
                           beta_init = beta_init0_all, weight = weight.count,
                           tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                           prop_presence = prop.presence, get_var = FALSE)
  
  rho_estimated_all = res.obs.wc.all$rho_estimates # the estimated correlation 

  # remove the OTU when executing Newton algorithm, since that the correlated freq.table.ref is 0.5 to avoid error in Newton algorithm
  cols_with_0_5 <- apply(freq.table.ref, 2, function(col) all(col == 0.5))
  cols_with_0_5_indices <- which(cols_with_0_5) # exactly the reference taxa
  columns_to_keep <- !apply(freq.table.ref, 2, function(col) all(col == 0.5))
  if (any(!columns_to_keep)) {
    freq.table.ref <- freq.table.ref[, columns_to_keep]
    N_denominator <- N_denominator[, columns_to_keep]
    beta_init0 <- beta_init0[, columns_to_keep]
    weight.count <- weight.count[, columns_to_keep]
    if(we) weight.equal<- weight.equal[, columns_to_keep]
  }
  dim(freq.table.ref)
  dim(N_denominator)
  dim(beta_init0)
  dim(weight.count)
  if(we) dim(weight.equal)
  
  if (!we.only) {
    res.obs.wc <- Newton_cor(freq_table = freq.table.ref, N_denominator = N_denominator, X = X,
                         beta_init = beta_init0, weight = weight.count,
                         tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                         prop_presence = prop.presence, get_var = FALSE)
    beta.wc <- matrix(NA, nrow = n.trait, ncol = n.taxa.common.filter)
    beta.wc[,-ref.taxon] <- res.obs.wc[[1]][1:n.trait,,drop=FALSE]
    beta.wc[,ref.taxon] <- 0
    #  which(is.nan(res.obs.wc$beta[1, ]))
  }
  
  if (we){
    res.obs.we <- Newton_cor(freq_table = freq.table.ref, N_denominator = N_denominator, X = X, #XX = XX
                         beta_init = beta_init0, weight = weight.equal,
                         tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                         prop_presence = prop.presence, get_var = FALSE)
    beta.we <- matrix(NA, nrow = n.trait, ncol = n.taxa.common.filter)
    beta.we[,-ref.taxon] <- res.obs.wc[[1]][1:n.trait,,drop=FALSE]
    beta.we[,ref.taxon] <- 0
  }
  
  if (n.taxa1.only.filter > 0) {
    res1.obs <- Newton_cor(freq_table = freq.table1.ref, N_denominator = N_denominator1, X = X1,
                       beta_init = array(0, dim = c(n.X1, n.taxa1.only.filter)), weight = weight1,
                       tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                       prop_presence = prop.presence1, get_var = FALSE)
    if (!we.only) beta.wc <- cbind(beta.wc, res1.obs[[1]][1:n.trait,,drop=FALSE] )
    beta.we <- cbind(beta.we, res1.obs[[1]][1:n.trait,,drop=FALSE] )
  }
  
  if (n.taxa2.only.filter > 0) {
    res2.obs <- Newton_cor(freq_table = freq.table2.ref, N_denominator = N_denominator2, X = X2, 
                       beta_init = array(0, dim = c(n.X2, n.taxa2.only.filter)), weight = weight2,
                       tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                       prop_presence = prop.presence2, get_var = FALSE)
    if (!we.only) beta.wc <- cbind(beta.wc, res2.obs[[1]][1:n.trait,,drop=FALSE] )
    beta.we <- cbind(beta.we, res2.obs[[1]][1:n.trait,,drop=FALSE])
  }
  
  # -----------------------
  # Remove the Taxa that did not converge
  # -----------------------
  
  Not_Converge_Taxa = which(res.obs.wc$final_step_sum < 3) # not include the reference taxon
  # consider the reference taxon's influence of the indices
  
  if (length(Not_Converge_Taxa) > 0) {
    Not_Converge_Taxa <- sapply(Not_Converge_Taxa, function(index) {
      if (index < ref.taxon) {
        return(index)
      } else {
        return(index + 1)
      }
    })
  } else {
    print("Not_Converge_Taxa is an empty vector (integer(0))")
  }
  
  Not_Converge_Taxa.df = as.data.frame(t(Not_Converge_Taxa))
  colnames(Not_Converge_Taxa.df) <- colnames(table.comTaxa.filter)[Not_Converge_Taxa]
  
  if (length(Not_Converge_Taxa) > 0) {
    if(clean == TRUE) { # remove the taxa that don't converge
    table.comTaxa.filter <- table.comTaxa.filter[,-Not_Converge_Taxa]
    
    n.sam.com <- nrow(table.comTaxa.filter)
    n.taxa.common.filter <- ncol(table.comTaxa.filter)
    
    # find the new reference OTU
    ref.taxon <- NULL
    if (is.null(ref.taxon)) {
      freq.table.comTaxa.filter <- table.comTaxa.filter/rowSums(table.comTaxa.filter)
      mean.freq <- colMeans(freq.table.comTaxa.filter)
      sd.freq <- colSds(freq.table.comTaxa.filter)
      ref.taxon <- which.max(mean.freq - sd.freq)
      rm(freq.table.comTaxa.filter)
    }
    table.comTaxa.filter[table.comTaxa.filter[, ref.taxon]==0, ref.taxon] <- 1
    
    freq.table.ref <- table.comTaxa.filter/(table.comTaxa.filter + table.comTaxa.filter[, ref.taxon]) # Y_j/(Y_j+Y_J)
    N_denominator = table.comTaxa.filter + table.comTaxa.filter[, ref.taxon] # Y_j+Y_J #weight.count
    dim(freq.table.ref)
    dim(N_denominator)
    
    if (!we.only) weight.count = (table.comTaxa.filter + table.comTaxa.filter[, ref.taxon])
    if(we) weight.equal = matrix(rep(colMeans(table.comTaxa.filter + table.comTaxa.filter[, ref.taxon]), n.sam.com), ncol=n.taxa.common.filter, byrow=TRUE) # temporary
    prop.presence <- colMeans(table.comTaxa.filter > 0) # rarity analysis
    beta_init0 <- array(0, dim = c(n.X, n.taxa.common.filter))
    
    cols_with_0_5 <- apply(freq.table.ref, 2, function(col) all(col == 0.5))
    cols_with_0_5_indices <- which(cols_with_0_5) # exactly the reference taxa
    columns_to_keep <- !apply(freq.table.ref, 2, function(col) all(col == 0.5))
    if (any(!columns_to_keep)) {
      freq.table.ref <- freq.table.ref[, columns_to_keep]
      N_denominator <- N_denominator[, columns_to_keep]
      beta_init0 <- beta_init0[, columns_to_keep]
      weight.count <- weight.count[, columns_to_keep]
      if(we) weight.equal<- weight.equal[, columns_to_keep]
    }
    dim(freq.table.ref)
    dim(N_denominator)
    dim(beta_init0)
    dim(weight.count)
    
    if (!we.only) {
      res.obs.wc <- Newton_cor(freq_table = freq.table.ref, N_denominator = N_denominator, X = X,
                           beta_init = beta_init0, weight = weight.count,
                           tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                           prop_presence = prop.presence, get_var = FALSE)
      beta.wc <- matrix(NA, nrow = n.trait, ncol = n.taxa.common.filter)
      beta.wc[,-ref.taxon] <- res.obs.wc[[1]][1:n.trait,,drop=FALSE]
      beta.wc[,ref.taxon] <- 0
    }
    
    if (we){
      res.obs.we <- Newton_cor(freq_table = freq.table.ref, N_denominator = N_denominator, X = X, #XX = XX
                           beta_init = beta_init0, weight = weight.equal,
                           tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                           prop_presence = prop.presence, get_var = FALSE)
      beta.we <- matrix(NA, nrow = n.trait, ncol = n.taxa.common.filter)
      beta.we[,-ref.taxon] <- res.obs.wc[[1]][1:n.trait,,drop=FALSE]
      beta.we[,ref.taxon] <- 0
    }
    
    # new results
    Not_Converge_Taxa = which(res.obs.wc$final_step_sum < 3) # don't include the reference taxon
    if (length(Not_Converge_Taxa) > 0) {
      Not_Converge_Taxa <- sapply(Not_Converge_Taxa, function(index) {
        if (index < ref.taxon) {
          return(index)
        } else {
          return(index + 1)
        }
      })
    } else {
      print("After cleaning, Not_Converge_Taxa is an empty vector (integer(0))")
    }
    Not_Converge_Taxa.df = as.data.frame(t(Not_Converge_Taxa))
    colnames(Not_Converge_Taxa.df) <- colnames(table.comTaxa.filter)[Not_Converge_Taxa]
  }
  }
  if (!we.only) beta.wc <- beta.wc - rowMedians(beta.wc)
  if(we) beta.we <- beta.we - rowMedians(beta.we)
  
  shrinkage <- ifelse(is.null(C), 0.5, 0.75)
  if (!we.only) beta_init_wc = rbind(matrix(0, nrow=n.trait, ncol=ncol(res.obs.wc[[1]])), shrinkage*res.obs.wc[[1]][(n.trait+1):n.X,,drop=FALSE])
  if(we) beta_init_we = rbind(matrix(0, nrow=n.trait, ncol=ncol(res.obs.we[[1]])), shrinkage*res.obs.we[[1]][(n.trait+1):n.X,,drop=FALSE])
  if (n.taxa1.only.filter > 0) beta_init1 = rbind(matrix(0, nrow=n.trait, ncol=ncol(res1.obs[[1]])), shrinkage*res1.obs[[1]][(n.trait+1):n.X1,,drop=FALSE])
  if (n.taxa2.only.filter > 0) beta_init2 = rbind(matrix(0, nrow=n.trait, ncol=ncol(res2.obs[[1]])), shrinkage*res2.obs[[1]][(n.trait+1):n.X2,,drop=FALSE])
  
  # -----------------------
  # Permutation
  # -----------------------
  
  if (is.null(seed)) seed = sample(1:10^6, 1)
  set.seed(seed)
  
  n.taxa <- n.taxa.common.filter + n.taxa1.only.filter + n.taxa2.only.filter # used to be n.taxa.union.filter
  
  if (is.null(n.perm.max)) n.perm.max = n.taxa * n.rej.stop * (1/fdr.nominal)
  
  if (!we.only) {
    beta.perm.wc <- array(NA, dim = c(n.trait, n.taxa, n.perm.max))
    n.rej.taxa.left.wc <- rep(0, n.taxa)
    n.rej.taxa.right.wc <- n.rej.taxa.left.wc
  }
  
  beta.perm.we <- array(NA, dim = c(n.trait, n.taxa, n.perm.max))
  n.rej.taxa.left.we <- rep(0, n.taxa)
  n.rej.taxa.right.we <- n.rej.taxa.left.we
  
  pnull.taxa.both <- array(NA, dim=c(n.trait*2, n.taxa, n.perm.max))
  
  tol.eq = 10^-8
  n.perm.block <- 1000
  n.block <- n.perm.max/n.perm.block
  n.perm.core <- n.perm.block/n.cores
  
  # Function Definition  ----------------------------------------------------
  
  if (!we.only) {
    parallel.perm.wc <- function(i) {
      res.perm.wc.fun <- perm_Newton_cor(freq_table = freq.table.ref, Yr = Yr, N_denominator = N_denominator, X = X,
                                     beta_init = beta_init_wc, weight = weight.count,
                                     perm = perm.mat[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                                     tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                                     prop_presence = prop.presence, get_var = FALSE)
      return(res.perm.wc.fun$beta)
    }
  }
  if (we){
    parallel.perm.we <- function(i) {
      res.perm.we.fun <- perm_Newton_cor(freq_table = freq.table.ref, Yr = Yr, N_denominator = N_denominator, X = X,
                                     beta_init = beta_init_we, weight = weight.equal,
                                     perm = perm.mat[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                                     tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                                     prop_presence = prop.presence, get_var = FALSE)
      return(res.perm.we.fun$beta)
    }
  }
  parallel.perm1 <- function(i) {
    res.perm1.tmp <- perm_Newton_cor(freq_table = freq.table1.ref, Yr = Y1r, N_denominator = N_denominator1, X = X1,
                beta_init = beta_init1, weight = weight1,
                perm = perm.mat1[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                prop_presence = prop.presence1, get_var = FALSE)
    return(res.perm1.tmp$beta)
  }
  parallel.perm2 <- function(i) {
    res.perm2.tmp <-perm_Newton_cor(freq_table = freq.table2.ref, Yr = Y2r, N_denominator = N_denominator2, X = X2,
                beta_init = beta_init2, weight = weight2,
                perm = perm.mat2[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                prop_presence = prop.presence2, get_var = FALSE)
    return(res.perm2.tmp$beta)
  }
  
  n.perm <- 0 # used to be n.perm.completed
  # used to save variance of each beta estimation across all slices in the first i.block
  variance_first_block = vector(mode="numeric")
  # For Loop to implement permutation  --------------------------------------
  for (i.block in 1:n.block) {
    
    # Permutation with three strata
    perm.mat.common <- t(shuffleSet(n.sam.common, n.perm.block))
    perm.mat <- rbind(perm.mat.common, perm.mat.common + n.sam.common)
    if (n.sam1.only > 0) {
      perm.mat1.only <- t(shuffleSet(n.sam1.only, n.perm.block))
      perm.mat <- rbind(perm.mat1.only, (perm.mat + n.sam1.only))
    }
    if (n.sam2.only > 0) {
      perm.mat2.only <- t(shuffleSet(n.sam2.only, n.perm.block))
      perm.mat <- rbind(perm.mat, perm.mat2.only + n.sam1.only + 2*n.sam.common)
    }
    perm.mat <- perm.mat - 1
    if (n.taxa1.only.filter > 0) {
      perm.mat1 <- perm.mat.common
      if (n.sam1.only > 0) {
        perm.mat1 <- rbind(perm.mat1.only, perm.mat.common + n.sam1.only)
      }
      perm.mat1 <- perm.mat1 - 1
    }
    if (n.taxa2.only.filter > 0) {
      perm.mat2 <- perm.mat.common
      if (n.sam2.only > 0) {
        perm.mat2 <- rbind(perm.mat.common, perm.mat2.only + n.sam.common)
      }
      perm.mat2 <- perm.mat2 -1
    }
    
    if (verbose) message(paste("permutations:", n.perm + 1))
    res.perm.wc = array(NA, dim = c(n.trait, n.taxa, dim(perm.mat)[2]))
    res.perm.we = array(NA, dim = c(n.trait, n.taxa, dim(perm.mat)[2]))
    
    if (n.cores > 1) {
      
      if (Sys.info()[['sysname']] == 'Windows') {
        if (!we.only) parallel.stat.wc = bplapply(0:(n.cores - 1), parallel.perm.wc, BPPARAM = MulticoreParam(workers=n.cores))
        if (we) parallel.stat.we = bplapply(0:(n.cores - 1), parallel.perm.we, BPPARAM = MulticoreParam(workers=n.cores))
        if (n.taxa1.only.filter > 0) parallel.stat.only.1 = bplapply(0:(n.cores - 1), parallel.perm1, BPPARAM = MulticoreParam(workers=n.cores))
        if (n.taxa2.only.filter > 0) parallel.stat.only.2 = bplapply(0:(n.cores - 1), parallel.perm2, BPPARAM = MulticoreParam(workers=n.cores))
      } else {
        if (!we.only) parallel.stat.wc = mclapply(0:(n.cores - 1), parallel.perm.wc, mc.cores = n.cores)
        if (we)  parallel.stat.we = mclapply(0:(n.cores - 1), parallel.perm.we, mc.cores = n.cores)
        if (n.taxa1.only.filter > 0) parallel.stat.only.1 = mclapply(0:(n.cores - 1), parallel.perm1, mc.cores = n.cores)
        if (n.taxa2.only.filter > 0) parallel.stat.only.2 = mclapply(0:(n.cores - 1), parallel.perm2, mc.cores = n.cores)
      }
      if (!we.only) {res.perm.wc.tmp <- abind(parallel.stat.wc)
      res.perm.wc[, -ref.taxon,] <- res.perm.wc.tmp
      res.perm.wc[, ref.taxon, ] <- 0 # beta estimation of the ref taxa
      }
      if (we) {res.perm.we.tmp <- abind(parallel.stat.we)
      res.perm.we[, -ref.taxon,] <- res.perm.we.tmp
      res.perm.we[, ref.taxon, ] <- 0 # beta estimation of the ref taxa
      }
      if (n.taxa1.only.filter > 0) res.perm1 <- abind(parallel.stat.only.1)
      if (n.taxa2.only.filter > 0) res.perm2 <- abind(parallel.stat.only.2)
      
    } else {
      if (!we.only) {
        res.perm.wc.tmp <- perm_Newton_cor(freq_table = freq.table.ref, Yr = Yr, N_denominator = N_denominator, X = X,
                                       beta_init = beta_init_wc, weight = weight.count,
                                       perm = perm.mat,
                                       tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                                       prop_presence = prop.presence, get_var = FALSE)
        res.perm.wc[, -ref.taxon,] <- res.perm.wc.tmp$beta
        res.perm.wc[, ref.taxon, ] <- 0 # beta estimation of the ref taxa
        
      }
      if (we) {
        res.perm.we.tmp <- perm_Newton_cor(freq_table = freq.table.ref, Yr = Yr,  N_denominator = N_denominator, X = X,
                                       beta_init = beta_init_we, weight = weight.equal,
                                       perm = perm.mat,
                                       tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = robust_var,
                                       prop_presence = prop.presence, get_var = FALSE)
        res.perm.we[, -ref.taxon,] <- res.perm.we.tmp$beta
        res.perm.we[, ref.taxon, ] <- 0 # beta estimation of the ref taxa
      }
      if (n.taxa1.only.filter > 0){
       res.perm1 <- perm_Newton_cor(freq_table = freq.table1.ref, Yr = Y1r, N_denominator = N_denominator1, X = X1,
                                                            beta_init = beta_init1, weight = weight1,
                                                            perm = perm.mat1,
                                                            tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                                            prop_presence = prop.presence1, get_var = FALSE)
      }
      if (n.taxa2.only.filter > 0) res.perm2 <- perm_Newton_cor(freq_table = freq.table2.ref, Yr = Y2r,  N_denominator = N_denominator2,  X = X2,
                                                            beta_init = beta_init2, weight = weight2,
                                                            perm = perm.mat2,
                                                            tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                                            prop_presence = prop.presence2, get_var = FALSE)
    }
    
    n.perm <- i.block*n.perm.block
    n.perm.inv <- 1 / (n.perm+1)
    
    # if (i.block ==1){
    #   for (i.var in 1:dim(res.perm.wc)[2]){
    #     variance_first_block[i.var] = var(res.perm.wc[,i.var,])
    #     #variance_per_block <- apply(res.perm.wc, 2, var)
    #   }
    # }
    # variance_df <- as.data.frame(t(variance_first_block))
    # colnames(variance_df) <- colnames(table.comTaxa.filter)
    
    if (!we.only) if (n.taxa1.only.filter > 0) res.perm.wc <- abind(res.perm.wc, res.perm1, along=2)
    if (!we.only) if (n.taxa2.only.filter > 0) res.perm.wc <- abind(res.perm.wc, res.perm2, along=2)
    if (we) if (n.taxa1.only.filter > 0) res.perm.we <- abind(res.perm.we, res.perm1, along=2)
    if (we) if (n.taxa2.only.filter > 0) res.perm.we <- abind(res.perm.we, res.perm2, along=2)
    
    w <- ((i.block-1)*n.perm.block+1):n.perm
    if (we) {
      for (i.trait in 1:n.trait) {
        beta.perm.we[i.trait,,w] <- t(t(res.perm.we[i.trait,,]) - colMedians(res.perm.we[i.trait,,], na.rm = TRUE))
      }
      
      diff <- sweep(beta.perm.we[,,w,drop=FALSE], c(1,2), beta.we, FUN="-")
      n.rej.taxa.equal.we <- 0.5*(abs(diff) < tol.eq) # 2 148 1000
      n.rej.taxa.left.we <- n.rej.taxa.left.we + rowSums( (diff < - tol.eq) + n.rej.taxa.equal.we, dims=2) # 2 148
      n.rej.taxa.right.we <- n.rej.taxa.right.we + rowSums( (diff > tol.eq) + n.rej.taxa.equal.we, dims=2) # 2 148
      n.rej.taxa.we <- 2*pmin(n.rej.taxa.left.we, n.rej.taxa.right.we) + 1 # 2 148
      if (n.trait >= 2) {
        pmin.taxa.we <- colMins(n.rej.taxa.we) # 148
        pnull.taxa.we <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.we[,,1:n.perm,drop=FALSE]), ties.method = "average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
        pnullmin.taxa.we <- array(colMins(pnull.taxa.we, dim.=c(n.trait, n.taxa*n.perm)), dim = c(n.taxa, n.perm)) # 148 1000
        diff <- pnullmin.taxa.we - c(pmin.taxa.we) # 148 1000
        n.rej.taxa.we2 <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
        p.taxa.we <- n.rej.taxa.we2 * n.perm.inv # 1 148
      } else {
        p.taxa.we <- n.rej.taxa.we * n.perm.inv # 1 148
      }
      q.taxa.we <- fdr.Sandve(p.taxa.we) # 148
    }
    
    if (!we.only) {
      
      for (i.trait in 1:n.trait) {
        beta.perm.wc[i.trait,,w] <- t(t(res.perm.wc[i.trait,,]) - colMedians(res.perm.wc[i.trait,,]))
      }
      diff <- sweep(beta.perm.wc[,,w,drop=FALSE], c(1,2), beta.wc, FUN="-")
      n.rej.taxa.equal.wc <- 0.5*(abs(diff) < tol.eq) # 2 148 1000
      n.rej.taxa.left.wc <- n.rej.taxa.left.wc + rowSums( (diff < - tol.eq) + n.rej.taxa.equal.wc, dims=2) # 2 148
      n.rej.taxa.right.wc <- n.rej.taxa.right.wc + rowSums( (diff > tol.eq) + n.rej.taxa.equal.wc, dims=2) # 2 148
      n.rej.taxa.wc <- 2*pmin(n.rej.taxa.left.wc, n.rej.taxa.right.wc) + 1 # 2 148
      if (n.trait >= 2) {
        pmin.taxa.wc <- colMins(n.rej.taxa.wc) # 148
        pnull.taxa.wc <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.wc[,,1:n.perm,drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
        pnullmin.taxa.wc <- array(colMins(pnull.taxa.wc, dim.=c(n.trait, n.taxa*n.perm)), dim = c(n.taxa, n.perm))   # 148 1000
        diff <- pnullmin.taxa.wc - c(pmin.taxa.wc)
        n.rej.taxa.wc2 <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
        p.taxa.wc <- n.rej.taxa.wc2 * n.perm.inv
      } else {
        p.taxa.wc <- n.rej.taxa.wc * n.perm.inv
      }
      q.taxa.wc <- fdr.Sandve(p.taxa.wc)
      
      if (we) pmin.taxa.omni <- colMins(rbind(n.rej.taxa.wc, n.rej.taxa.we)) # 148
      if (n.trait == 1) {
        if (we) pnull.taxa.we <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.we[,,1:n.perm,drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
        pnull.taxa.wc <- n.perm + 0.5 - array(rowRanks(abs(beta.perm.wc[,,1:n.perm,drop=FALSE]), ties.method="average", dim.=c(n.trait*n.taxa, n.perm)), dim = c(n.trait, n.taxa, n.perm)) # 2 148 1000
      }
      pnull.taxa.both[1:n.trait,,1:n.perm] <- pnull.taxa.wc
      if (we) {pnull.taxa.both[(n.trait+1):(2*n.trait),,1:n.perm] <- pnull.taxa.we # 4 148 1000
      pnullmin.taxa.omni <- array(colMins(pnull.taxa.both[,,1:n.perm], dim.=c(n.trait*2, n.taxa*n.perm)), dim = c(n.taxa, n.perm)) # 148 1000
      diff <- pnullmin.taxa.omni- c(pmin.taxa.omni)
      n.rej.taxa.omni <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
      p.taxa.omni <- n.rej.taxa.omni * n.perm.inv
      q.taxa.omni <- fdr.Sandve(p.taxa.omni)
      }
    }
    
    if (!we.only) {
      if (n.perm >= 10000
          & all(q.taxa.wc <= fdr.nominal | n.rej.taxa.wc >= n.rej.stop)
          # & all(q.taxa.we <= fdr.nominal | n.rej.taxa.we >= n.rej.stop) #  if (we)
          # & all(q.taxa.omni <= fdr.nominal | n.rej.taxa.omni >= n.rej.stop)
      ) break
    } else {
      if (n.perm >= 10000
          & all(q.taxa.we <= fdr.nominal | n.rej.taxa.we >= n.rej.stop)
      ) break
    }
    
  } # permutation
  
  taxa.union <- colnames(table.comTaxa.filter)
  if (n.taxa1.only.filter > 0) taxa.union <- c(taxa.union, colnames(table1.only.filter))
  if (n.taxa2.only.filter > 0) taxa.union <- c(taxa.union, colnames(table2.only.filter))
  
  # Multiple test correction by BH
  if (we) {q.taxa.we <- p.adjust(p.taxa.we, method = "BH")
  detected.taxa.we <- taxa.union[which(q.taxa.we < fdr.nominal)]}
  if (!we.only) {
    q.taxa.wc <- p.adjust(p.taxa.wc, method = "BH")
    # q.taxa.omni <- p.adjust(p.taxa.omni, method = "BH")
    detected.taxa.wc <- taxa.union[which(q.taxa.wc < fdr.nominal)]
    # detected.taxa.omni <- taxa.union[which(q.taxa.omni < fdr.nominal)]
  }
  
  # ------------------------
  # Global p-value
  # ------------------------
  if (we){
    if (n.trait == 1) {
      beta.all.we <- cbind(beta.we[1,], beta.perm.we[1,,1:n.perm])
    } else {
      beta.all.we <- cbind(pmin.taxa.we, pnullmin.taxa.we)  # 148 1001
    }
    
    r <- rowRanks(beta.all.we, ties.method = "average") # 148 1001
    p.taxa1.we <- (2*pmin(r, n.perm+2-r) - 1.5) # omit n.perm.inv # 148 1001
    stat.global.we <- colSums(1/p.taxa1.we)
    p.global.we <- (sum(stat.global.we[1] <= stat.global.we[-1]) + 1) * n.perm.inv
  }
  if (!we.only) {
    if (n.trait == 1) {
      beta.all.wc <- cbind(beta.wc[1,], beta.perm.wc[1,,1:n.perm])
    } else {
      beta.all.wc <- cbind(pmin.taxa.wc, pnullmin.taxa.wc)
    }
    r <- rowRanks(beta.all.wc, ties.method = "average") # 148 1001
    p.taxa1.wc <- (2*pmin(r, n.perm+2-r) - 1.5)
    stat.global.wc <- colSums(1/p.taxa1.wc)
    p.global.wc <- (sum(stat.global.wc[1] <= stat.global.wc[-1]) + 1) * n.perm.inv
    if (we){
      beta.all.omni <- cbind(pmin.taxa.omni, pnullmin.taxa.omni) # 148 1001
      p.taxa1.omni <- (rowRanks(beta.all.omni, ties.method = "average") - 0.5) # 148 1001
      stat.global.omni <- colSums(1/p.taxa1.omni) # 1001
      p.global.omni<- (sum(stat.global.omni[1] <= stat.global.omni[-1]) + 1) * n.perm.inv
    }
  }
  if (we){
    beta.we <- matrix(beta.we, nrow=n.trait)
    colnames(beta.we) <- taxa.union
    p.taxa.we <- matrix(c(p.taxa.we), nrow = 1)
    q.taxa.we <- matrix(c(q.taxa.we), nrow = 1)
    colnames(p.taxa.we) <- taxa.union
    colnames(q.taxa.we) <- taxa.union
  }
  
  if (!we.only) {
    beta.wc <- matrix(beta.wc, nrow=n.trait)
    colnames(beta.wc) <- taxa.union
    p.taxa.wc <- matrix(c(p.taxa.wc), nrow = 1)
    q.taxa.wc <- matrix(c(q.taxa.wc), nrow = 1)
    colnames(p.taxa.wc) <- taxa.union
    colnames(q.taxa.wc) <- taxa.union
    
    # p.taxa.omni <- matrix(c(p.taxa.omni), nrow = 1) #  if (we)
    # q.taxa.omni <- matrix(c(q.taxa.omni), nrow = 1)
    # colnames(p.taxa.omni) <- taxa.union
    # colnames(q.taxa.omni) <- taxa.union
  } else {
    beta.wc <- NULL
    p.taxa.wc <- NULL
    q.taxa.wc <- NULL
    detected.taxa.wc <- NULL
    p.global.wc <- NULL
    
    p.taxa.omni <- NULL
    q.taxa.omni <- NULL
    detected.taxa.omni = NULL
    p.global.omni = NULL
  }
  
  return(list(effect.size.wc = beta.wc,
              p.taxa.wc = p.taxa.wc,
              q.taxa.wc = q.taxa.wc,
              detected.taxa.wc = detected.taxa.wc,
              p.global.wc = p.global.wc,
              Not_Converge_Taxa.df = Not_Converge_Taxa.df,
              freq.table.ref = freq.table.ref,
              rho_estimated_all = rho_estimated_all,
              
              # variance_df = variance_df,
              # effect.size.we = beta.we,
              # p.taxa.we = p.taxa.we,
              # q.taxa.we = q.taxa.we,
              # detected.taxa.we = detected.taxa.we,
              # p.global.we = p.global.we,
              
              # p.taxa.omni = p.taxa.omni,
              # q.taxa.omni = q.taxa.omni,
              # detected.taxa.omni = detected.taxa.omni,
              # p.global.omni = p.global.omni,
              
              n.perm.completed = n.perm,
              ref.taxon = ref.taxon,
              seed = seed))
} # Com2seq


fdr.Sandve = function(p.otu) {
    m = length(p.otu)
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)
    qval.orig = qval.sort[mat]
    results = qval.orig
    
    return(results)
    
} # fdr.Sandve