#' Soft Thresholding
#'
#' This internal function conduct soft-thresholding on a square matrix
#' @param x Matrix that is thresholded
#' @param t Soft-thresholding criterion
#' @param shrink_diag Bool indicator on whether to shrink the diagonal terms of X.
#' @return thresholded x
#' @noRd
soft_thresh=function(x,t,shrink_diag=FALSE){
  if(!shrink_diag){
    y=sign(x)*ifelse(abs(x)-t>0,abs(x)-t,0)
    diag(y)=diag(x)
    return(y)
  }else{
    return(sign(x)*ifelse(abs(x)-t>0,abs(x)-t,0))
  }
}
