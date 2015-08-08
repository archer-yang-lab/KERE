plot.KER <- function(x, 
    color = FALSE, label = FALSE, ...) {
    alpha <- x$alpha[-1,]
    lambda <- x$lambda
    index <- log(lambda)
    iname <- "Log Lambda"
    xlab <- iname
    ylab <- "Coefficients"
    dotlist <- list(...)
    type <- dotlist$type
    if (is.null(type)) {
        if (color == FALSE) 
            matplot(index, t(alpha), lty = 1, xlab = xlab, ylab = ylab, 
                type = "l", pch = 500, col = gray.colors(12, 
                  start = 0.05, end = 0.7, gamma = 2.2), ...) else matplot(index, t(alpha), lty = 1, xlab = xlab, ylab = ylab, 
            type = "l", pch = 500, ...)
    } else matplot(index, t(alpha), lty = 1, xlab = xlab, ylab = ylab, 
        ...)
    if (label) {
        nnz <- length(which)
        xpos <- max(index)
        pos <- 4
        if (xvar == "lambda") {
            xpos <- min(index)
            pos <- 2
        }
        xpos <- rep(xpos, nnz)
        ypos <- alpha[, ncol(alpha)]
        text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }
} 
