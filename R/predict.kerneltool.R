predict.kerneltool <- function(object, newx, type = c("link"), ...) {
    type <- match.arg(type)
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
    nfit
} 
