source("kernels.R", chdir = TRUE)
source("source.R") # R code for generating the FHT data

nobs = 100
p = 200
rho = 0.9

X=genx2(nobs,p,rho)
y=genjerry(X,3)
## initialize kernel function
rbf <- polydot(degree = 1, scale = 1, offset = 1)

## calculate kernel matrix
Kmat <- kernelMatrix(rbf,X)


ulam  <-  c(0.4,0.5,0.6,0.7)
omega  <-  0.3
omega <- as.double(omega)
mbd <- 2 * max(1-omega, omega)
nlam <- length(ulam)
eps <- 10e-14
maxit  <-  1e6

dl <- function(r, omega)
{
	d = 2.0*ifelse(r>0, omega, 1-omega)*r	
}


B <- alpmat
for (l in 1:nlam)
{
	ri <- y-(cbind(1,Kmat)%*%B[,l])
	L = dl(ri,omega)
	yxl <- - rbind(1,Kmat) %*% L / nobs + 2 * ulam[l] * cbind(0,rbind(0,Kmat)) %*% alpmat[, l]
	if (crossprod(yxl)>1e-6) print(paste("this is",crossprod(yxl)," at ",l))
}
