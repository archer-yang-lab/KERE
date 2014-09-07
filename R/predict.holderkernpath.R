predict.holderkernpath2 <- function(object, kern, x, newx,
                                 type = c("class", "link"), ...) {
  type <- match.arg(type)
  nfit <- kernelMult(kern, newx, x, object$alpha[-1, ])
  nfit <- sweep(nfit, MARGIN = 2, object$alpha[1, ], "+")
  switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
} 