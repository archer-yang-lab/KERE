err <- function(n, maxit) {
    if (n == 0) 
        msg <- ""
    if (n < 0) {    
        msg <- paste("convergence for ", -n, "th lambda value not reached after maxit=", 
            maxit, " iterations; solutions for larger lambdas returned", 
            sep = "")
        n <- -1
        msg <- paste("From kerneltool fortran code -", msg)
    }
    list(n = n, msg = msg)
}