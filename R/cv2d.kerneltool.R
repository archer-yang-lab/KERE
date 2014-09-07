cv2d.kerneltool <- function(x, y, kname = "rbfdot", lambda = NULL, sigma = NULL, ...) {
		if (is.null(sigma)) {
	        stop("user must provide a sigma sequence")}
		nsigma <- length(sigma)
		cputime <- rep(NA, nsigma)
		for(i in seq.int(nsigma)){
			kern <- do.call(kname, list(sigma[i]))
			cputime[i] <- system.time(cv_out <- cv.kerneltool(x, y, kern,					     
					lambda = lambda, ...))[3]
			which <- (cvout$lambda == cvout$lambda.min)				
			}
} 
