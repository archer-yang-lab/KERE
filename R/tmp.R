source("kernels.R", chdir = TRUE)
X = matrix(rnorm(10*10,0,1),10,10)

## initialize kernel function
rbf <- rbfdot(sigma = 0.05)
rbf
## calculate kernel matrix
kernelMatrix(rbf, X)

Kmat <- kernelMatrix(rbf,X)
