#' Self-Concordant Proximal Newton Method(primal) for Graphical Lasso
#'
#' This function uses local Hessian estimate and primal solution to search for direction and step size in proximal gradient method and solve for inverse covariance matrix
#' @param S Sample covariance matrix
#' @param n Sample size
#' @param lambda Penalty parameter
#' @param warm.start Bool indicator for whether to do warm start for estimation
#' @param init.Omega Initial value for inverse covariance matrix, only needed when warm.start=TRUE
#' @return Estimated inverse covariance matrix
#' @export
Newton_Prox_primal_glasso=function(S,n,lambda,warm.start=F,init.Omega=NULL){
  p=nrow(S)
  sigma=(5-sqrt(17))/4
  if(warm.start&&(!is.null(init.Omega))){
    Omega0=Omega1=init.Omega
  }else{
    Omega0=Omega1=diag(p)
  }
  thresh=1e-5
  lam=1
  iter=0
  prim=diag(p)

  while (lam>thresh) {
    iter=iter+1
    Linv=solve(chol(Omega0))
    Omega0.inv=Linv%*%t(Linv)
    prim=solve_primal(Omega0,Omega0.inv,S,prim,lambda)
    delta=prim-Omega0
    lam=sqrt(sum(diag(Omega0.inv%*%delta%*%Omega0.inv%*%delta)))
    alpha=ifelse(lam>sigma,1/(1+lam),1)
    Omega1=Omega0+alpha*delta
    Omega0=Omega1

  }
  return(Omega0)
}
