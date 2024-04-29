#' Fully Symmetric Inverse
#'
#' This internal function derive the inverse of a matrix that is numerically symmetric
#' @param mat Matrix that is to be inverted
#' @return inverse of mat
#' @noRd
cool_inverse=function(mat){
  inv_mat=solve(mat)
  inv_mat=(inv_mat+t(inv_mat))/2
  return(inv_mat)
}
