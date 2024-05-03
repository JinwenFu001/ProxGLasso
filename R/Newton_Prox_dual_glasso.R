#' Self-Concordant Proximal Newton Method(dual) for Graphical Lasso
#'
#' This function uses local Hessian estimate and dual solution to search for direction and step size in proximal gradient method and solve for inverse covariance matrix
#' @param S Sample covariance matrix
#' @param n Sample size
#' @param lambda Penalty parameter
#' @param warm.start Bool indicator for whether to do warm start for estimation
#' @param init.Omega Initial value for inverse covariance matrix, only needed when warm.start=TRUE
#' @return Estimated inverse covariance matrix
#' @examples
#' n=10
#' p=5
#' set.seed(42)
#' samples=Generate_AR1_pair(n,p,rho=0.7,cov_AR1 = TRUE)
#' S=cov(samples$X)/n*(n-1)
#' SC_Newton_dual_est=Newton_Prox_dual_glasso(S,n,lambda = 0.1)
#' @export
Newton_Prox_dual_glasso=function(S,n,lambda,warm.start=FALSE,init.Omega=NULL){
  p=nrow(S)
  sigma=(5-sqrt(17))/4
  if(warm.start&&(!is.null(init.Omega))){
    Omega0=Omega1=init.Omega
  }else{
    Omega0=Omega1=diag(p)
  }
  U=diag(p)
  thresh=1e-10
  lam=1
  iter=0

  while (lam>thresh) {
    iter=iter+1
    Q=(Omega0%*%S%*%Omega0-2*Omega0)/lambda
    Q=force_symmetric(Q)
    U=solve_dual(Omega0,Q,U)
    W=Omega0%*%(S+lambda*U)
    lam=sqrt(abs(p-2*sum(diag(W))+sum(diag(W%*%W))))
    delta=-(Omega0%*%(S+lambda*U)%*%t(Omega0)-Omega0)
    delta=force_symmetric(delta)
    alpha=ifelse(lam>sigma,1/(1+lam),1)
    Omega1=Omega0+alpha*delta
    Omega0=Omega1*(abs(Omega1)>1e-5)
  }
  return(Omega0)
}
