kerneltool <- function(x, y, nlambda = 100, method = c( 
    "logit", "ls"), lambda.factor = 0.01, 
    lambda = NULL, exclude, standardize = TRUE, 
    eps = 1e-08, maxit = 1e+06, delta = 2, omega = 0.5) {
    #################################################################################
    #data setup
    method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    if (length(y) != nobs) 
        stop("x and y have different number of observations")
    #################################################################################
    #parameter setup
    maxit <- as.integer(maxit)
    isd <- as.integer(standardize)
    eps <- as.double(eps)
    if (!missing(exclude)) {
        jd <- match(exclude, seq(nvars), 0)
        if (!all(jd > 0)) 
            stop("Some excluded variables out of range")
        jd <- as.integer(c(length(jd), jd))
    } else jd <- as.integer(0)
    #################################################################################
    #lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
        if (lambda.factor >= 1) 
            stop("lambda.factor should be less than 1")
        flmin <- as.double(lambda.factor)
        ulam <- double(1)
    } else {
        #flmin=1 if user define lambda
        flmin <- as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    #################################################################################
    fit <- switch(method, 
	hhsvm = hsvmpath(x, y, nlam, flmin, 
        ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, 
        nobs, nvars, vnames), 
	er = erpath(x, y, nlam, flmin, 
        ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, omega, 
        nobs, nvars, vnames))
    if (is.null(lambda)) 
        fit$lambda <- lamfix(fit$lambda)
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("kerneltool", class(fit))
    fit
} 
