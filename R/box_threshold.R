#' Box Threshold
#'
#' This internal function threshold X inside a box
#' @param X Matrix or vector or scalar we're interested in
#' @param bound The bound that specifies the box
#' @return Thresholded X
#' @noRd
box_threshold=function(X,bound=1){
  (abs(X)-bound<= 0)*(X)+(X< -bound)*(-bound)+(X>bound)*bound
}
