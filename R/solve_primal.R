#' Solve Primal Problem with FISTA
#'
#' This internal function solves the primal subproblem of proximal Newton method
#' @param Omega inverse covariance in last step
#' @param Omega.inv inverse of last estimate of inverse covariance matrix
#' @param S sample covariance matrix
#' @param prim solution in last step
#' @param lambda penalty parameter
#' @return new primal solution
#' @noRd
solve_primal=function(Omega,Omega.inv, S,prim,lambda){
  p=nrow(Omega)
  Y0=Y1=X0=X1=prim
  lam=1/eigen(Omega.inv)$values[1]^2
  t0=t1=1

  f.deriv=function(U){
    return(Omega.inv%*%U%*%Omega.inv+S-2*Omega.inv)
  }


  obj=function(U){
    sum((S-2*Omega.inv)*U)+sum((Omega.inv%*%U%*%Omega.inv%*%U))/2+lambda*sum(abs(U))
  }

  iter=0
  thresh=1e-5
  dis=1
  while(dis>thresh){
    iter=iter+1
    temp_y=Y0-lam*f.deriv(Y0)
    X1=soft_thresh(temp_y,lambda*lam)
    t1=(1+sqrt(1+4*t0^2))/2
    Y1=X1+(t0-1)/t1*(X1-X0)
    dis=sqrt(sum((Y1-Y0)^2))/p^2

    Y0=Y1
    X0=X1
    t0=t1
  }
  return(X1)
}
