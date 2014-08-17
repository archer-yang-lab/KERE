## S4 object definitions and assigment/accessor functions for the slots.
##
## created  10.09.03 alexandros karatzoglou
## updated  23.08.05 


setClass("kernel",representation("function",kpar="list"))
setClass("kernelMatrix",representation("matrix"),prototype=structure(.Data=matrix()))

setGeneric("as.kernelMatrix",function(x, center = FALSE) standardGeneric("as.kernelMatrix"))
setMethod("as.kernelMatrix", signature(x = "matrix"),
function(x, center = FALSE)
{

  if(center){
    m <- dim(x)[1]
    x <- t(t(x - colSums(x)/m) -  rowSums(x)/m) + sum(x)/m^2
  }
  
  return(new("kernelMatrix",.Data = x))
})




fun <- function(object) standardGeneric("kpar")
setGeneric("kpar",fun)
setMethod("kpar","kernel", function(object) object@kpar)





## kernel functions
## Functions for computing a kernel value, matrix, matrix-vector
## product and quadratic form
##
## author : alexandros karatzoglou


## Define the kernel objects,
## functions with an additional slot for the kernel parameter list.
## kernel functions take two vector arguments and return a scalar (dot product)


rbfdot<- function(sigma=1)
  {

    rval <- function(x,y=NULL)
    {
       if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y))))  
        # sigma/2 or sigma ??
      }
    }
     return(new("rbfkernel",.Data=rval,kpar=list(sigma=sigma)))
  }
setClass("rbfkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

laplacedot<- function(sigma=1)
  {
    rval <- function(x,y=NULL)
    {
       if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(-sigma*sqrt(-(round(2*crossprod(x,y) - crossprod(x) - crossprod(y),9)))))
      }
    }
     return(new("laplacekernel",.Data=rval,kpar=list(sigma=sigma)))
  }

setClass("laplacekernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

besseldot<- function(sigma = 1, order = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        lim <- 1/(gamma(order+1)*2^(order))
        bkt <- sigma*sqrt(-(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
        if(bkt < 10e-5)
          res <- lim
        else
          res <- besselJ(bkt,order)*(bkt^(-order))
        return((res/lim)^degree)
      }
    }
     return(new("besselkernel",.Data=rval,kpar=list(sigma=sigma ,order = order ,degree = degree)))
  }

setClass("besselkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

anovadot<- function(sigma = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")

          res <- sum(exp(- sigma * (x - y)^2))
        return((res)^degree)
      }
    }
     return(new("anovakernel",.Data=rval,kpar=list(sigma=sigma ,degree = degree)))
  }

setClass("anovakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))


splinedot<- function()
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        minv <- pmin(x,y)
        res <- 1 + x*y*(1+minv) - ((x+y)/2)*minv^2 + (minv^3)/3
          fres <- prod(res)
        return(fres)
      }
    }
     return(new("splinekernel",.Data=rval,kpar=list()))
  }

setClass("splinekernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))



fourierdot <- function(sigma = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
	res <- 	(1 - sigma^2)/2*(1 - 2*sigma*cos(x - y) + sigma^2)
	fres <- prod(res)
        return(fres)
      }
    }
     return(new("fourierkernel",.Data=rval,kpar=list()))
  }

setClass("fourierkernel",prototype=structure(.Data=function(){},kpar=list(sigma = 1)),contains=c("kernel"))





tanhdot <- function(scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a  vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        tanh(scale*crossprod(x)+offset)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        tanh(scale*crossprod(x,y)+offset)
      }
    }
  return(new("tanhkernel",.Data=rval,kpar=list(scale=scale,offset=offset)))
}
setClass("tanhkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

setClass("polykernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

polydot <- function(degree = 1, scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        (scale*crossprod(x)+offset)^degree
        }

      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        (scale*crossprod(x,y)+offset)^degree
          }

      }
  return(new("polykernel",.Data=rval,kpar=list(degree=degree,scale=scale,offset=offset)))
}

setClass("vanillakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

vanilladot <- function( )
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        crossprod(x)
        }

      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        crossprod(x,y)
          }

      }
  return(new("vanillakernel",.Data=rval,kpar=list()))
}



## show method for kernel functions

setMethod("show",signature(object="kernel"),
          function(object)
          {
            switch(class(object),
                   "rbfkernel" = cat(paste("Gaussian Radial Basis kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
		   "laplacekernel" = cat(paste("Laplace kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
                   "besselkernel" = cat(paste("Bessel kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"order = ",kpar(object)$order, "degree = ", kpar(object)$degree,"\n")),
                    "anovakernel" = cat(paste("Anova RBF kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma, "degree = ", kpar(object)$degree,"\n")),
                   "tanhkernel" = cat(paste("Hyperbolic Tangent kernel function.", "\n","Hyperparameters :","scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "polykernel" = cat(paste("Polynomial kernel function.", "\n","Hyperparameters :","degree = ",kpar(object)$degree," scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "vanillakernel" = cat(paste("Linear (vanilla) kernel function.", "\n")),
                   "splinekernel" = cat(paste("Spline kernel function.", "\n")), 
                  )
          })







## Functions that return usefull kernel calculations (kernel matrix etc.)

## kernelMatrix function takes two or three arguments

kernelMatrix <- function(kernel, x, y=NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(x,"matrix")) stop("x must be a matrix")
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  n <- nrow(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  if(is.null(y)){
    for(i in 1:n) {
      for(j in i:n) {
        res1[i,j]  <- kernel(x[i,],x[j,])
      }
    }
    res1 <- res1 + t(res1)
    diag(res1) <- diag(res1)/2
    
  
 }
  if (is(y,"matrix")){
    m<-dim(y)[1]
    res1 <- matrix(0,dim(x)[1],dim(y)[1])
    for(i in 1:n) {
      for(j in 1:m) {
        res1[i,j] <- kernel(x[i,],y[j,])
      }
    }
  }

  return(as.kernelMatrix(res1))
}

setGeneric("kernelMatrix",function(kernel, x, y = NULL) standardGeneric("kernelMatrix"))



kernelMatrix.rbfkernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(2*sigma*(res[i,] - dota - rep(dota[i],n)))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="rbfkernel"),kernelMatrix.rbfkernel)

kernelMatrix.laplacekernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(-sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9)))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(-sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9)))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="laplacekernel"),kernelMatrix.laplacekernel)

kernelMatrix.besselkernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  nu = kpar(kernel)$order
  ni = kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n){
      xx <- sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9))
      res[i,] <- besselJ(xx,nu)*(xx^(-nu))
      res[i,which(xx<10e-5)] <- lim
    }
    return(as.kernelMatrix((res/lim)^ni))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m){
      xx <- sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9))
      res[,i] <- besselJ(xx,nu)*(xx^(-nu))
      res[which(xx<10e-5),i] <- lim
    }     
    return(as.kernelMatrix((res/lim)^ni))
  }
}
setMethod("kernelMatrix",signature(kernel="besselkernel"),kernelMatrix.besselkernel)


kernelMatrix.anovakernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  degree = kpar(kernel)$degree
  n <- dim(x)[1]
  if (is(x,"matrix") && is.null(y)){
    a <- matrix(0,  dim(x)[2], n)
    res <- matrix(0, n ,n)
    for (i in 1:n)
      {
        a[rep(TRUE,dim(x)[2]), rep(TRUE,n)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(a - t(x))^2))^degree
      }
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    
    m <- dim(y)[1]
    b <- matrix(0, dim(x)[2],m)
    res <- matrix(0, dim(x)[1],m)
    for( i in 1:n)
      {
        b[rep(TRUE,dim(x)[2]), rep(TRUE,m)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(b - t(y))^2))^degree
      }
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="anovakernel"),kernelMatrix.anovakernel)


kernelMatrix.polykernel <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  scale = kpar(kernel)$scale
  offset = kpar(kernel)$offset
  degree = kpar(kernel)$degree
  if (is(x,"matrix") && is.null(y))
    {
      res <- (scale*crossprod(t(x))+offset)^degree
      return(as.kernelMatrix(res))
    }
      if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- (scale*crossprod(t(x),t(y)) + offset)^degree
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="polykernel"),kernelMatrix.polykernel)

kernelMatrix.vanilla <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- crossprod(t(x),t(y))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="vanillakernel"),kernelMatrix.vanilla)

kernelMatrix.tanhkernel <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if (is(x,"matrix") && is.null(y)){
    scale = kpar(kernel)$scale
    offset = kpar(kernel)$offset
    res <- tanh(scale*crossprod(t(x)) + offset)
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- tanh(scale*crossprod(t(x),t(y)) + offset)
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="tanhkernel"),kernelMatrix.tanhkernel)


kernelMatrix.splinekernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  degree = kpar(kernel)$degree
  n <- dim(x)[1]
  if (is(x,"matrix") && is.null(y)){
    a <- matrix(0,  dim(x)[2], n)
    res <- matrix(0, n ,n)
	x <- t(x)
    for (i in 1:n)
      {
	dr <- 	x + x[,i]	
	dp <-   x * x[,i]
	dm <-   pmin(x,x[,i])
	res[i,] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
      }
        return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    
    m <- dim(y)[1]
    b <- matrix(0, dim(x)[2],m)
    res <- matrix(0, dim(x)[1],m)
	x <- t(x)
	y <- t(y)
    for( i in 1:n)
      {
	dr <- 	y + x[,i]	
	dp <-   y * x[,i]
	dm <-   pmin(y,x[,i])
	res[i,] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
   	}
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="splinekernel"),kernelMatrix.splinekernel)

