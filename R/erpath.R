erpath <- function(x, y, nlam, flmin, ulam, isd, 
    eps, jd, maxit, omega, nobs, nvars, vnames) {
    #################################################################################
    #data setup
    y <- as.double(y)
    if (omega <= 0 || omega >= 1) 
        stop("omega must be in (0,1)")
	omega <- as.double(omega)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("erlassoNET", omega, nobs, nvars, as.double(x), 
        as.double(y), jd, nlam, flmin, ulam, 
        eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
        beta = double(nobs * nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), 
        PACKAGE = "kerneltool")
    #################################################################################
    # output
    outlist <- fit
    class(outlist) <- c("erpath")
    outlist
} 
