cv2d.KERE <-
  function(x, y, kname = "rbfdot", lambda = NULL, sigma = NULL, ...) {
    if (is.null(sigma)) {
      stop("user must provide a sigma sequence")
    }
    mm.cvm <- Inf
    for (i in seq.int(length(sigma))) {
      kern <- do.call(kname, list(sigma[i]))
      cv_out <- cv.KERE(x, y, kern, lambda = lambda, ...)
      if (mm.cvm > cv_out$cvm.min) {
        mm.cvm <- cv_out$cvm.min
        mm.lambda <- cv_out$lambda.min
        loc.sigma <- i
      }
      cat("sigma ", i, " completed.\n")
    }
    loc.lambda <- which(mm.lambda == lambda)
    list(
      mm.cvm = mm.cvm, loc.lambda = loc.lambda,
      loc.sigma = loc.sigma, mm.lambda = mm.lambda,
      mm.sigma = sigma[loc.sigma]
    )
  }
