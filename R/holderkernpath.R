holderkernpath <- function(x, y, Kmat, nlam, ulam, 
    eps, maxit, qval, nobs) {
    #################################################################################
    #data setup
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    if (!all(y %in% c(-1, 1))) 
        stop("y should be a factor with two levels")
    if (qval < 0) 
        stop("qval must be positive")
    qval <- as.double(qval)
    eigen_result <- eigen(Kmat, symmetric = TRUE)
    Umat <- eigen_result$vectors
    Dvec <- eigen_result$values
    Ksum <- colSums(Kmat)
    Bmat <- - Ksum %o% Ksum / nobs
    ################################################################################
   	#call Fortran core
    fit <- .Fortran("holderkern", qval, 
			as.double(Kmat), as.double(Umat),
			as.double(Dvec), as.double(Ksum), as.double(Bmat), 
			nobs, as.double(y), nlam, ulam, eps, maxit, anlam = integer(1), 
			npass = integer(nlam), jerr = integer(1), 
			alpmat = double((nobs+1) * nlam),
			PACKAGE = "kerneltool")
    ################################################################################
    # output
    anlam <- fit$anlam
    alpha <- matrix(fit$alpmat[seq((nobs+1) * anlam)], nobs+1, anlam) 
    outlist <- list(alpha = alpha, lambda = ulam[seq(anlam)], npass = fit$npass, jerr = fit$jerr)
    class(outlist) <- c("holderkernpath")
    outlist
} 
