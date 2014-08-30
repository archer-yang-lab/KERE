expdirect <- function(x, y, Kmat, nlam, ulam, 
    eps, maxit, omega, nobs) {
    #################################################################################
    #data setup
    y <- as.double(y)
    if (omega <= 0 || omega >= 1) 
        stop("omega must be in (0,1)")
	omega <- as.double(omega)
	
	Ksum <- colSums(Kmat)
	#################################################################################
	K0 = rbind(0, cbind(0, Kmat))
    Ki = rbind(1, Kmat)
	mbd <- 2 * max(1-omega, omega)
	npass <- rep(0,nlam)
    r <- y # r = 0 in classification case
	alpmat <- matrix(0, nobs+1, nlam)
	alpvec = rep(0, nobs+1)
	for(l in 1:nlam) {
		dif <- rep(NA, nobs+1)
		# computing Ku inverse
		Amat <- Kmat %*% Kmat + 2*nobs*ulam[l]*Kmat/mbd
		KU <- matrix(NA, nobs+1, nobs+1)
		KU[1, 1] <- nobs
		KU[1, 2:(nobs+1)]  <- Ksum
		KU[2:(nobs+1), 1]  <- Ksum
		KU[2:(nobs+1), 2:(nobs+1)]  <- Amat
		KUinv <- solve(KU)
		# for debug
		# KUtmp = rbind(c(nobs, Ksum), cbind(Ksum, Kmat%*%Kmat + 2*nobs*ulam[l]*Kmat / mbd))
		# KUinv1  <- solve(KUtmp)
		# update alpha
		oldalpvec = alpvec
		while(1){
		    phi <- 2.0*ifelse(r>0, omega, 1-omega)*r
		    oldalpvec <- alpvec
		    alpvec <-  oldalpvec + (2*nobs/mbd) * KUinv %*% (-ulam[l]*K0%*%oldalpvec + 0.5*Ki%*%phi/nobs) 
		    alpvec <-  drop(alpvec)
		    dif <- alpvec - oldalpvec
		    r <- r - dif %*% Ki
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
    outlist <- list(alpha = alpmat)
    class(outlist) <- c("expdirect")
    outlist
} 
