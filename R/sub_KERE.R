sub.KERE <- function(x, y, kern, lambda = NULL,
                    nsub = 10, ...) {
  #####################################
  #data setup
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  nobs <- NROW(x)
  #Subsampling
  effective_size <- nobs-nobs%%nsub
  sub_size <- effective_size/nsub
  fitmat <- array(NA, c(sub_size+1, length(lambda), nsub))
  subid <- sample(rep(seq(nsub), length = effective_size))
  for (i in seq(nsub)) {
    which <- subid == i
    y_sub <- y[which]
    fit <- KERE(
      x = x[which, , drop = FALSE],
      y = y_sub, kern = kern, lambda = lambda, omega = omega, ...
    )
	fitmat[,,i] <- fit$alpha
  }
  alpha <- apply(fitmat,c(1,2),mean)
  alpha
}
