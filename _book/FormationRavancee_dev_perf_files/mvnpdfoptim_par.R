#' @rdname mvnpdf
#' @importFrom itertools isplitIndices
#' @importFrom future plan multisession availableCores
#' @importFrom future.apply future_sapply
#' @export
#' @examples
#' \dontrun{
#' n <- 10000
#' mb <- microbenchmark::microbenchmark(
#'   mvtnorm::dmvnorm(matrix(1.96, nrow = n, ncol = 2)),
#'   mvnpdfsmart(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE),
#'   mvnpdfoptim_par(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE),
#'   times=10L)
#'mb
#'}
#'
mvnpdfoptim_par <- function(x, mean =  rep(0, nrow(x)), varcovM = diag(nrow(x)),
                            Log=TRUE){

  if(!is.matrix(x)){
    x <- matrix(x, ncol=1)
  }

  n <- ncol(x)
  p <- nrow(x)
  x0 <- x-mean

  Rinv = backsolve(chol(varcovM), x=diag(p))

  y_par <- future.apply::future_sapply(X = 1:p, FUN = function(i) {
    xRinv <- apply(X=x0[, i, drop=FALSE], MARGIN=2, FUN=crossprod, y=Rinv)
    logSqrtDetvarcovM <- sum(log(diag(Rinv)))

    quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
    y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)

    if(!Log){
      y <- exp(y)
    }
    return(y)
  })
  res <- list(x = x, y = y_par)
  class(res) <- "mvnpdf"
  return(res)
}


#' @export
mvnpdfoptim_parIter <- function(x, mean =  rep(0, nrow(x)),
                                varcovM = diag(nrow(x)), Log=TRUE, ncores = 1){

  if(!is.matrix(x)){
    x <- matrix(x, ncol=1)
  }

  n <- ncol(x)
  p <- nrow(x)
  x0 <- x-mean

  Rinv = backsolve(chol(varcovM), x=diag(p))
  iter <- itertools::isplitIndices(n = p, chunks = ncores)
  y_par <- future.apply::future_sapply(X = iter, FUN = function(i) {
    xRinv <- apply(X=x0[, i, drop=FALSE], MARGIN=2, FUN=crossprod, y=Rinv)
    logSqrtDetvarcovM <- sum(log(diag(Rinv)))

    quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
    y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)

    if(!Log){
      y <- exp(y)
    }
    return(y)
  })
  res <- list(x = x, y = y_par)
  class(res) <- "mvnpdf"
  return(res)
}
