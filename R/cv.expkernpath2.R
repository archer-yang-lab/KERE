cv.expkernpath2 <- function(outlist, lambda, x, y, kern, foldid, 
    pred.loss, qval, omega) {
    typenames <- c(misclass = "Misclassification Error", loss = "Expectile Loss")
    if (pred.loss == "default") 
        pred.loss <- "loss"
    if (!match(pred.loss, c("loss"), FALSE)) {
        warning("Only 'loss' available for expectile regression; 'loss' used")
        pred.loss <- "loss"
    }

    y <- as.double(y)
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, kern, x[!which, ,drop = FALSE], x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- ercls(y-predmat, omega)
	cvob <- cvcompute(cvraw, foldid, nlams)
	cvraw <- cvob$cvraw
	N <- cvob$N
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
