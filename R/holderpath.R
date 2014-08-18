hsvmpath <- function(x, y, nlam, flmin, ulam, 
    eps, jd, maxit, delta, nobs) {
    #################################################################################
    #data setup
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    if (!all(y %in% c(-1, 1))) 
        stop("y should be a factor with two levels")
    if (delta < 0) 
        stop("delta must be non-negative")
    delta <- as.double(delta)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("hsvmlassoNET", delta, nobs, nvars, 
        as.double(x), as.double(y), jd, nlam, 
        flmin, ulam, eps, maxit, nalam = integer(1), b0 = double(nlam), 
        beta = double(nobs * nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), 
        PACKAGE = "kerneltool")
    #################################################################################
    # output
    outlist <- fit
    class(outlist) <- c("hsvmpath")
    outlist
} 
