#' Is Positive Definite
#'
#' This internal function tells whether or not a matrix is positive definite
#' @param A Matrix interested in
#' @return Bool value of whether or not A is positive definite
#' @noRd
is.pd <- function(A) {
  # Use tryCatch to attempt the Cholesky decomposition
  tryCatch({
    chol(A)
    TRUE  # Return TRUE if chol succeeds
  }, error = function(e) {
    FALSE  # Return FALSE if an error occurs
  })
}
