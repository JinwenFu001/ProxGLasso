#' Simulate data points with AR1 covariance structure
#'
#' This function simulate data points follow multivariate normal distribution with given dimension and structure
#' @importFrom SimDesign rmvnorm
#' @param n Sample size
#' @param p Sample dimension
#' @param rho Coefficient of AR1 model
#' @param cov_AR1 Bool indicator for whether the covariance or inverse covariance has AR1 structure
#' @return A list containing true covariance, inverse covariance and simulated data
#' @examples
#' n=10
#' p=5
#' set.seed(42)
#' samples=Generate_AR1_pair(n,p,rho=0.7,cov_AR1 = TRUE)
#' @export
Generate_AR1_pair=function(n,p,rho,cov_AR1=TRUE){
  mat1=matrix(0,nrow = p,ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      mat1[i,j]=rho^abs(i-j)
    }
  }
  Linv=solve(chol(mat1))
  mat2=Linv%*%t(Linv)
  mat2=ifelse(abs(mat2)>1e-10,mat2,0)
  if(cov_AR1){
    cov_mat=mat1
    prec_mat=mat2
  }else{
    cov_mat=mat2
    prec_mat=mat1
  }
  X=SimDesign::rmvnorm(n,sigma=cov_mat)
  return(list(covariance=cov_mat,precision=prec_mat,X=X))
}
