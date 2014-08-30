err <- function(n, maxit) {
    if (n == 0) 
        msg <- ""
    if (n < 0) {    
        msg <- paste("convergence for ", -n, "th lambda value not reached after maxit=", 
            maxit, " iterations; solutions for larger lambdas returned", 
            sep = "")
        n <- -1
        msg <- paste("From kerneltool fortran code -", msg)
    }
    list(n = n, msg = msg)
}


cvcompute <- function(mat, foldid, nlams) {
    ###Computes the weighted mean and SD within folds, and
    #   hence
    #   the se of the mean
    nfolds <- max(foldid)
    outmat <- matrix(NA, nfolds, ncol(mat))
    good <- matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] <- NA
    for (i in seq(nfolds)) {
        mati <- mat[foldid == i, ]
        outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
        good[i, seq(nlams[i])] <- 1
    }
    N <- apply(good, 2, sum)
    list(cvraw = outmat, N = N)
}


hubercls <- function(r, delta) {
  (1 - r - 0.5 * delta) * (r <= 1 - delta) + 0.5 * (1 - r)^2 / delta * (r <= 1) * (r > 1 - delta)
} 


ercls <- function(r, omega) {
  abs(omega - (r < 0)) * r^2
}


getmin <- function(lambda, cvm, cvsd) {
    cvmin <- min(cvm)
    idmin <- cvm <= cvmin
    lambda.min <- max(lambda[idmin])
    idmin <- match(lambda.min, lambda)
    semin <- (cvm + cvsd)[idmin]
    idmin <- cvm <= semin
    lambda.1se <- max(lambda[idmin])
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

error.bars <- function(x, upper, lower, width = 0.02, 
    ...) {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}






