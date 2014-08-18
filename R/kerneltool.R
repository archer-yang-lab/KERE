kerneltool <- function(x, y, Kmat, method = c("hhsvm", "er", "exp","holder"), 					     
	lambda = NULL, eps = 1e-08, maxit = 1e+06, delta = 2, 
	omega = 0.5, gamma = 1e-06) {
    #################################################################################
    #data setup
    method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
	Kmat <- as.matrix(Kmat)
	diag(Kmat) <- diag(Kmat) + gamma 
    np <- dim(x)
    nobs <- as.integer(np[1])
    if (length(y) != nobs) 
        stop("x and y have different number of observations")
    #################################################################################
    #parameter setup
    maxit <- as.integer(maxit)
    eps <- as.double(eps)
    #################################################################################
    #lambda setup
    if (is.null(lambda)) {
        stop("user must provide a lambda sequence")
    } else {
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    #################################################################################
    fit <- switch(method, 
	holder = holderkernpath(x, y, Kmat, nlam, ulam, eps, maxit, delta, 
        nobs), 
	exp = expkernpath(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
        nobs),
	hhsvm = hsvmpath(x, y, Kmat, nlam, ulam, eps, maxit, delta, 
	        nobs),
	er = erpath(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
        nobs))
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("kerneltool", class(fit))
    fit
} 
