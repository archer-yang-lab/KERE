hsvmpath <- function(x, y, Kmat, nlam, ulam, 
    eps, maxit, delta, nobs) {
    #################################################################################
    #data setup
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    if (!all(y %in% c(-1, 1))) 
        stop("y should be a factor with two levels.")
    if (delta <= 0) 
        stop("delta must be positive.")
	delta <- as.double(delta)
  eigen_result <- eigen(Kmat, symmetric = TRUE)
	Umat <- eigen_result$vectors
	Dvec <- eigen_result$values
	Ksum <- colSums(Kmat)
	Bmat <- - Ksum %o% Ksum / nobs
	#################################################################################
	K0 = rbind(0, cbind(0, Kmat))
  Ki = rbind(1, Kmat)
	mbd = (delta+1.0)*(delta+1.0)/delta
  minv = 1 / mbd
	decib = delta / (delta + 1.0)
	fdr = - decib ** (delta + 1.0)
	npass <- rep(0,nlam) # for each lambda, compute the number of iterates
  r <- rep(0,nobs) # r_i = y_i (\beta_0 + K_i \alpha)
	alpmat <- matrix(0, nobs+1, nlam) # save alpha values for all values lambdas
	alpvec = rep(0, nobs+1)
	for(l in 1:nlam) {
		dif <- rep(NA, nobs+1)
		# computing Ku inverse
		Ainv <- Umat %*% diag(1/(Dvec^2 + 2*nobs*ulam[l]*Dvec*minv)) %*% t(Umat)
		BAmat <- Bmat %*% Ainv
		Ginv <- 1 / (1 + sum(diag(BAmat)))
		Qinv <- Ainv - Ginv * Ainv %*% BAmat
		QKsum <- Qinv %*% Ksum / nobs
		Mtmp <- (1 + crossprod(QKsum,Ksum)) / nobs
		KUinv <- matrix(NA, nobs+1, nobs+1)
		KUinv[1, 1] <- Mtmp
		KUinv[1, 2:(nobs+1)]  <- -QKsum
		KUinv[2:(nobs+1), 1]  <- -QKsum
		KUinv[2:(nobs+1), 2:(nobs+1)]  <- Qinv
		# for debug
		# KUtmp = rbind(c(nobs, Ksum), cbind(Ksum, Kmat%*%Kmat + 2*nobs*ulam[l]*Kmat / mbd))
		# KUinv1  <- solve(KUtmp)
		# update alpha
		oldalpvec = alpvec
		while(1){
		    phi <- ifelse(r>decib, r**(-delta-1)*fdr, -1.0)
		    oldalpvec <- alpvec
		    alpvec <-  oldalpvec - (nobs/mbd) * KUinv %*% (2*ulam[l]*K0%*%oldalpvec + Ki%*%(y*phi)/nobs) 
		    alpvec <-  drop(alpvec)
		    dif <- alpvec - oldalpvec
		    r <- r + y * dif %*% Ki
		    r <- drop(r)
		    if(sum(dif^2)/sum(oldalpvec^2) < eps) break
		    npass[l] = npass[l] + 1
		    if(sum(npass) > maxit) break
		}
		alpmat[, l] <- alpvec
		if(sum(npass) > maxit) {
			break
			jerr = -l
		}
	}
    ################################################################################
    # output
    outlist <- list(alpha = alpmat, npass=npass)
    class(outlist) <- c("hsvmpath")
    outlist
} 
