#' Force Symmetric
#'
#' This internal function forces a square matrix to be numerically symmetric
#' @param X Matrix interested in
#' @return symmetric version of X
#' @noRd
force_symmetric=function(X){
  (X+t(X))/2
}
