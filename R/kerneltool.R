kerneltool <- function(x, y, kern, method = c("exp", "holder"), 					     
	lambda = NULL, eps = 1e-08, maxit = 1e+06, qval = 2, 
	omega = 0.5, gamma = 1e-06) {
    #################################################################################
    #data setup
    method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
	  Kmat <- kernelMatrix(kern,x)
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
	# expdirect = expdirect(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
	#         nobs),
	# expfast = expfast(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
	#         nobs),
    # exp = expkernpath(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
    #                   nobs),
    exp = expkernpath(x, y, Kmat, nlam, ulam, eps, maxit, omega, 
                      nobs),
    holder = holderkernpath(x, y, Kmat, nlam, ulam, eps, maxit, qval, 
                            nobs) 
    #holderexp = holderexppath(x, y, Kmat, nlam, ulam, eps, maxit, qval, 
    #                 nobs)
    )
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("kerneltool", class(fit))
    fit
} 
