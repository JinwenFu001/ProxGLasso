#' Backtracking Proximal Gradient Method for Graphical Lasso
#'
#' This function uses backtracking to search for step size in proximal gradient method and solve for inverse covariance matrix
#' @param S Sample covariance matrix
#' @param n Sample size
#' @param lambda Penalty parameter
#' @param warm.start Bool indicator for whether to do warm start for estimation
#' @param init.Omega Initial value for inverse covariance matrix, only needed when warm.start=TRUE
#' @return Estimated inverse covariance matrix
#' @export
BTProx_glasso=function(S,n,lambda,warm.start=F,init.Omega=NULL){
  p=nrow(S)
  S.eigen.max=eigen(S)$values[1]
  if(warm.start&&(!is.null(init.Omega))){
    Omega0=Omega1=init.Omega
  }else{
    Omega0=Omega1=diag(p)
  }
  mu0=mu1=1
  rate=0.5
  dis=1
  thresh=1e-7
  iter=0
  max.iter=1e5
  objval=1

  f=function(Omega){
    return(sum(S*Omega)-log(det(Omega)))
  }

  derivative.f=function(Omega.inv){
    return(S-Omega.inv)
  }

  obj_eval=function(S,Omega,lambda){
    return(sum(S*Omega)-log(det(Omega))+lambda*sum(abs(Omega-diag(diag(Omega)))))
  }

  Q_f=function(Omega,Omegat,Omegat.deriv,mut){
    return(f(Omegat)-f(Omega)+sum((Omega-Omegat)*Omegat.deriv)+sum((Omega-Omegat)^2)/2/mut)
  }


  while ((dis>thresh)&&(iter<=max.iter)) {
    iter=iter+1
    Linv=solve(chol(Omega0))
    Omega0.inv=Linv%*%t(Linv)
    mu.safe=min(min(diag(Linv))^2,n/(p*lambda+n*S.eigen.max))/n
    mu1=5
    Omega1=soft_thresh(Omega0-mu1*derivative.f(Omega0.inv),mu1*lambda)
    iter0=0
    while (( (!is.pd(Omega1))||(Q_f(Omega1,Omega0,derivative.f(Omega0.inv),mu1)<0))&iter0<=10) {
      iter0=iter0+1
      mu1=mu1*rate
      Omega1=soft_thresh(Omega0-mu1*derivative.f(Omega0.inv),mu1*lambda)
    }
    if(iter0==11){
      Omega1=soft_thresh(Omega0-mu.safe*derivative.f(Omega0.inv),mu.safe*lambda)
    }
    dis=abs(objval-obj_eval(S,Omega1,lambda))/abs(objval)
    objval=obj_eval(S,Omega1,lambda)
    Omega0=Omega1
  }
  return(Omega1)
}
