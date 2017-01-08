cv2d_H.KERE <-
  function(x, y, kname = "tanhdot", lambda = NULL, scale = NULL, ...) {
    if (is.null(scale)) {
      stop("user must provide a scale sequence")
    }
    mm.cvm <- Inf
    for (i in seq.int(length(scale))) {
      kern <- do.call(kname, list(scale[i], offset = 1))
      cv_out <- cv.KERE(x, y, kern, lambda = lambda, ...)
      if (mm.cvm > cv_out$cvm.min) {
        mm.cvm <- cv_out$cvm.min
        mm.lambda <- cv_out$lambda.min
        loc.scale <- i
      }
      cat("scale ", i, " completed.\n")
    }
    loc.lambda <- which(mm.lambda == lambda)
    list(
      mm.cvm = mm.cvm, loc.lambda = loc.lambda,
      loc.scale = loc.scale, mm.lambda = mm.lambda,
      mm.scale = scale[loc.scale]
    )
  }
