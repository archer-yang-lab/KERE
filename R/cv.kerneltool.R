cv.kerneltool <- function(x, y, kern, lambda = NULL, 
	pred.loss = c("loss", "misclass"), nfolds = 5, foldid, qval = 2.0, omega = 0.5, ...) {
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    x <- as.matrix(x)
    N <- NROW(x)
    # predict -> coef
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- kerneltool(x = x[!which, , drop = FALSE], 
            y = y_sub, kern = kern, lambda = lambda, 
			qval = qval, omega = omega, ...)
		cat("fold ", i, " completed.\n")
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(outlist[[1]])[[2]], sep = ".")
    lambda <- outlist[[1]]$lambda
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, kern, foldid, 
        pred.loss, qval, omega))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlo = cvm - cvsd, name = cvname)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.kerneltool"
    obj
} 
