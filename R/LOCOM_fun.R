#' A logistic regression model for testing differential abundance in compositional microbiome data (LOCOM)
#'
#' This function allows you to test
#' (1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on the log ratio of relative abundances between pairs of taxa, and
#' (2). whehter the whole community is associated with the trait (a global test).
#' The tests accommodate continuous, discrete (binary, categorical), and multivariate traits, and allows adjustment of confounders.
#'
#' This function uses a sequential stopping criterion (similar to that of Sandve et al. 2011) for the permutation procedure,
#' which stops when all taxon-level tests have either reached the pre-specified
#' number (default 100) of rejections or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) that is below the
#' nominal FDR level (default 0.2). The permutation procedure is always terminated if a pre-specified maximum number (see description of \code{n.perm.max}) of
#' permutations have been generated. The global test is based on all permutation replicates when the procedure stops/terminates.
#'
#'
#' @param otu.table the OTU table (or taxa count table) in which the rows correspond to samples and the columns correspond to OTUs (taxa).
#' @param Y the trait of interest. It can be an array, matrix, or data frame; a factor should be represented by the corresponding design matrix. All components in a matrix or data frame will be tested for microbial association jointly. 
#' @param C the other (confounding) covariates to be adjusted. See other requirements in \code{Y}.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param seed a user-supplied integer seed for the random number generator in the
#'   permutation procedure. The default is NULL, which means that an integer seed will be
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case the user wants to reproduce the permutation replicates.
#' @param n.perm.max the maximum number of permutations. The default is NULL, in which case \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal})
#'   is used where \code{n.otu} is the total number of OTUs (that have non-zero counts in at least one sample).
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation test
#'   statistic exceeds the observed test statistic) to obtain before stopping the permutation procedure. The
#'   default is 100.
#' @param filter.thresh a real value between 0 and 1; OTUs present in fewer than \code{filter.thresh} * 100\% samples are filtered out. The default is 0.2.
#' @param n.cores the number of cores to be used for parallel computing. The default is 1.
#' @param permute a logical value indicating whether to perform permutation. The default is TRUE. 
#' @return A list consisting of
#' \itemize{
#'   \item effect.size - effect size at each OTU, i.e., beta_j,1 - median_j'=1,...J(beta_j',1)
#'   \item p.otu - p-values for OTU-specific tests
#'   \item q.otu - q-values (adjusted p-values by the BH procedure) for OTU-specific tests
#'   \item detected.otu - detected OTUs (using the column names of the OTU table) at the nominal FDR based on \code{q.otu}
#'   \item p.global - p-value for the global test
#'   \item n.perm.completed - number of permutations completed
#'   \item seed - the seed used to generate the permutation replicates
#'   \item p.otu.asymptotic - asymptotic p-values for OTU-specific tests
#'   \item q.otu.asymptotic - asymptotic q-values (adjusted p-values by the BH procedure) for OTU-specific tests
#'   \item detected.otu.asymptotic - detected OTUs (using the column names of the OTU table) at the nominal FDR based on \code{q.otu.asymptotic}
#' }
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom parallel mclapply
#' @importFrom permute shuffleSet
#' @importFrom stats cor pchisq p.adjust lm resid
#' @importFrom matrixStats rowMedians colMedians colMins colSds rowRanks
#' @importFrom abind abind
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib LOCOM
#' @export
#' @examples
#' library(LOCOM)
#' data("throat.otu.table")
#' data("throat.meta")
#' data("throat.otu.taxonomy")
#' 
#' # filtering out three samples with antibiotic use
#' filter.out.sam <- which(throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")
#' throat.otu.table.filter <- throat.otu.table[-filter.out.sam,]
#' throat.meta.filter <- throat.meta[-filter.out.sam,]
#' Y <- ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
#' C <- ifelse(throat.meta.filter$Sex == "Male", 0, 1)
#' 
#' # filtering out rare OTUs
#' prop.presence <- colMeans(throat.otu.table.filter > 0)
#' filter.out.otu <- which(prop.presence < 0.4) # In practice, use threshold 0.2; 
#'                                              # here, filter out more OTUs for speed
#' if (length(filter.out.otu) > 0) {
#'     throat.otu.table.filter <- throat.otu.table.filter[, -filter.out.otu]
#'     prop.presence <- prop.presence[-filter.out.otu]
#' }
#' 
#' # running locom
#' res <- locom(otu.table = throat.otu.table.filter, Y = Y, C = C, 
#'              seed = 123, n.cores = 1, n.perm.max = 1000)
#' res$p.global 
#' res$detected.otu


locom <- function(otu.table, Y, C = NULL,
                  fdr.nominal = 0.2, seed = NULL,
                  n.perm.max = NULL, n.rej.stop = 100, n.cores = 4,
                  filter.thresh = 0.2, permute = TRUE){
    
    
    Firth.thresh = 0.4
    tol = 1e-6
    iter_max = 50
    
    if (is.data.frame(otu.table)) otu.table <- data.matrix(otu.table)
    if (is.data.frame(Y)) Y <- data.matrix(Y)
    if (is.data.frame(C)) C <- data.matrix(C)
    
    n.trait <- length(Y)/nrow(otu.table)
    Y <- matrix(Y, ncol=n.trait)
    if (!is.null(C)) {
        n.conf <- length(C)/nrow(otu.table)
        C <- matrix(C, ncol=n.conf)
    }
    
    # remove rare OTUs
    
    n.sam <- nrow(otu.table)
    w = which(colSums(otu.table>0)>= ceiling(filter.thresh * n.sam))
    if (length(w) < ncol(otu.table)) {
        cat(paste(ncol(otu.table)-length(w), ' OTU(s) present in fewer than ', ceiling(filter.thresh * n.sam), ' samples are filtered out', sep=""), "\n")
        otu.table = otu.table[,w,drop=FALSE]
    }
    n.otu <- ncol(otu.table)
    
    # find reference OTU
    
    ref.otu <- NULL
    if (is.null(ref.otu)) {
        mean.freq <- colMeans(otu.table/rowSums(otu.table))
        ref.otu <- which.max(mean.freq)
    }
    otu.table[which(otu.table[,ref.otu]==0), ref.otu] <- 1 # replace 0 by 1
    
    # -----------------------
    # Observed statistic
    # -----------------------
    
    if (is.null(C) == TRUE) {
        X <- cbind(Y, 1)
        Yr <- Y
    } else {
        X <- cbind(Y, 1, C)
        Yr <- resid(lm(Y ~ C))
        Yr <- matrix(Yr, ncol=length(Yr)/nrow(Y))
    }
    
    freq.table.ref <- otu.table/(otu.table + otu.table[,ref.otu])
    n.X <- ncol(X)
    XX <- CalculateXX(X)
    weight = (otu.table + otu.table[,ref.otu])
    prop.presence <- colMeans(otu.table>0)
    beta_init0 <- array(0, dim = c(n.X, n.otu))

    res.obs <- tryCatch(
        {
            Newton(freq_table = freq.table.ref, X = X, XX = XX, 
                   beta_init = beta_init0, weight = weight,
                   tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                   prop_presence = prop.presence, get_var = TRUE)
        },
        error=function(cond) {
            message(cond)
            return(NA)
        }
    )
    if (any(is.na(res.obs)) == TRUE) {
        return(NA)
    }    
    
    beta <- res.obs[[1]][1:n.trait,,drop=FALSE] 
    var.median <- apply(beta, 1, mckean.schrader.var)
    beta <- beta - rowMedians(beta)
    p.otu.asymptotic <- rep(NA, n.otu)
    for (j in 1:n.otu) {
        if (n.trait >=2) {
            p.otu.asymptotic[j] <- t(beta[,j,drop=FALSE]) %*% solve(res.obs[[2]][1:n.trait,1:n.trait,j] + diag(var.median)) %*% beta[,j,drop=FALSE]
        } else {
            p.otu.asymptotic[j] <- t(beta[,j,drop=FALSE]) %*% solve(res.obs[[2]][1:n.trait,1:n.trait,j] + var.median) %*% beta[,j,drop=FALSE]
        }
    }
    p.otu.asymptotic <- pchisq(p.otu.asymptotic, df=n.trait, lower.tail=FALSE)
    q.otu.asymptotic <- p.adjust(p.otu.asymptotic, method = "BH")
    detected.otu.asymptotic <- colnames(otu.table)[which(q.otu.asymptotic < fdr.nominal)]
    
    # -----------------------
    # Permutation
    # -----------------------
    
    p.otu <- NULL
    q.otu <- NULL
    
    if (permute) {
        
        shrinkage <- ifelse(is.null(C), 0.5, 0.75)
        beta_init = rbind(matrix(0, nrow=n.trait, ncol=n.otu), shrinkage*res.obs[[1]][(n.trait+1):n.X,,drop=FALSE])
    
        if (is.null(seed)) seed = sample(1:10^6, 1)
        set.seed(seed)
        
        if (is.null(n.perm.max)) n.perm.max = n.otu * n.rej.stop * (1/fdr.nominal)
        
        beta.perm <- array(NA, dim = c(n.trait, n.otu, n.perm.max))
        n.rej.otu.left <- 0
        n.rej.otu.right <- 0
        
        tol.eq = 10^-8
        n.perm.block <- 1000
        n.block <- n.perm.max/n.perm.block
        n.perm.core <- n.perm.block/n.cores
        
        
        parallel.perm <- function(i) {
            tryCatch(
                {
                    perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX, 
                                beta_init = beta_init, weight = weight,
                                perm = perm.mat[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                                tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                prop_presence = prop.presence, get_var = FALSE)
                },
                error=function(cond) {
                    return(list(res = NA, error =cond))
                }
            )    
        } # parallel.perm
    
        
        n.perm <- 0
        
        for (i.block in 1:n.block){
            
            perm.mat <- t(shuffleSet(n.sam, n.perm.block)) - 1
            
            cat("permutations:", n.perm + 1, "\n")
            
            if (n.cores > 1) {
                
                if (Sys.info()[['sysname']] == 'Windows') {
                    parallel.stat = bplapply(0:(n.cores - 1), parallel.perm, BPPARAM = MulticoreParam(workers=n.cores))
                } else {
                    parallel.stat = mclapply(0:(n.cores - 1), parallel.perm, mc.cores = n.cores)
                }
                
                # check whether there is any error in a single core
                error.list <- unlist(lapply(parallel.stat, function(x) any(is.na(x) == TRUE)))
                if (sum(error.list) > 0) {
                    error.idx <- min(which(error.list == TRUE))
                    message(parallel.stat[[error.idx]][2])
                    return(NA)
                }
                res.perm <- do.call(abind, parallel.stat) # default is to bind along the last dimension
                
            } else {
                res.perm <- tryCatch(
                    {
                        perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX,
                                    beta_init = beta_init, weight = weight,
                                    perm = perm.mat,
                                    tol = tol, iter_max = iter_max, Firth_thresh = Firth.thresh, robust_var = FALSE,
                                    prop_presence = prop.presence, get_var = FALSE)
                    },
                    error=function(cond) {
                        message(cond)
                        return(NA)
                    }
                )
                if (any(is.na(res.perm)) == TRUE) {
                    return(NA)
                }
            }
            
            n.perm <- i.block*n.perm.block
            n.perm.inv <- 1 / (n.perm+1)
            
            w <- ((i.block-1)*n.perm.block+1):n.perm
            for (i.trait in 1:n.trait) {
                beta.perm[i.trait,,w] <- t(t(res.perm[i.trait,,]) - colMedians(res.perm[i.trait,,]))
            }
            diff <- sweep(beta.perm[,,w,drop=FALSE], c(1,2), beta, FUN="-")
            n.rej.otu.equal <- 0.5*(abs(diff) < tol.eq) # 2 148 1000
            n.rej.otu.left <- n.rej.otu.left + rowSums( (diff < - tol.eq) + n.rej.otu.equal, dims=2) # 2 148
            n.rej.otu.right <- n.rej.otu.right + rowSums( (diff > tol.eq) + n.rej.otu.equal, dims=2) # 2 148
            n.rej.otu <- 2*pmin(n.rej.otu.left, n.rej.otu.right) + 1 # 2 148
            
            if (n.trait >= 2) {
                pmin.otu <- colMins(n.rej.otu) # 148
                pnull.otu <- n.perm + 0.5 - array(rowRanks(abs(beta.perm[,,1:n.perm,drop=FALSE]), ties.method = "average", dim.=c(n.trait*n.otu, n.perm)), dim = c(n.trait, n.otu, n.perm)) # 2 148 1000
                pnullmin.otu <- array(colMins(pnull.otu, dim.=c(n.trait, n.otu*n.perm)), dim = c(n.otu, n.perm)) # 148 1000 
                diff <- pnullmin.otu - c(pmin.otu) # 148 1000 
                n.rej.otu <- rowSums( (diff < - tol.eq) + 0.5 * ( abs(diff) < tol.eq) )
            }
            p.otu <- n.rej.otu * n.perm.inv   
            q.otu <- fdr.Sandve(p.otu)
            
            if (all(q.otu <= fdr.nominal | n.rej.otu >= n.rej.stop)) break
            
        } # permutation
        
        q.otu <- p.adjust(p.otu, method = "BH")
        detected.otu <- colnames(otu.table)[which(q.otu < fdr.nominal)]
        
        # ------------------------
        # Global p-value
        # ------------------------
        
        if (n.trait == 1) {
            beta.all <- cbind(beta[1,], beta.perm[1,,1:n.perm])
        } else {
            beta.all <- cbind(pmin.otu, pnullmin.otu)
        }
        r <- rowRanks(beta.all, ties.method = "average")
        p.otu1 <- (2*pmin(r, n.perm+2-r) - 1.5) # omit n.perm.inv # 148 1001
        
        stat.global <- colSums(1/p.otu1)
        p.global <- (sum(stat.global[1] <= stat.global[-1]) + 1) * n.perm.inv
    
    } # if (permute)
    
    otu.names <- colnames(otu.table)
    
    beta <- matrix(beta, nrow = n.trait)
    colnames(beta) <- otu.names
    p.otu.asymptotic <- matrix(p.otu.asymptotic, nrow = 1)
    q.otu.asymptotic <- matrix(q.otu.asymptotic, nrow = 1)
    colnames(p.otu.asymptotic) <- otu.names
    colnames(q.otu.asymptotic) <- otu.names

    if (!is.null(p.otu)) {
        p.otu <- matrix(p.otu, nrow = 1)
        q.otu <- matrix(q.otu, nrow = 1)
        colnames(p.otu) <- otu.names
        colnames(q.otu) <- otu.names
    }
    
    return(list(effect.size = beta,
                p.otu = p.otu,
                q.otu = q.otu,
                detected.otu = detected.otu,
                p.global = p.global,
                n.perm.completed = n.perm,
                seed = seed,
                p.otu.asymptotic = p.otu.asymptotic,
                q.otu.asymptotic = q.otu.asymptotic,
                detected.otu.asymptotic = detected.otu.asymptotic))
    
} # locom


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


mckean.schrader.var = function( x ) {
    
    # The mckean-schrader standard error and variance of a sample median
    
    n.data = length(x)
    c = max(1, round( (n.data+1)/2 - 0.98*sqrt(n.data) ) )    #  0.98=1.96/sqrt(4)
    x.sort = sort(x)
    med.se = ( x.sort[n.data-c+1] - x.sort[c] )/3.92     #    3.92=2*1.96
    med.var = med.se^2
    
    return(med.var)
    
} # mckean.schrader.var

