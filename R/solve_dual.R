#' Solve Dual Problem with FISTA
#'
#' This internal function solves the dual subproblem of proximal Newton method
#' @param Omega inverse covariance in last step
#' @param Q parameter
#' @param U solution in last step
#' @return new dual solution
#' @noRd
solve_dual=function(Omega,Q,U){
  p=nrow(Omega)
  Y0=Y1=X0=X1=U
  Omega_sqr=Omega%*%Omega
  lam=1/(eigen(Omega)$values[1])^2
  t0=t1=1

  f.deriv=function(U){
    return(Omega%*%U%*%Omega+Q)
  }

  obj=function(U){
    sum(Q*U)+sum(Omega_sqr*(U%*%U))/2
  }

  iter=0
  thresh=1e-5
  dis=1
  while(dis>thresh){
    iter=iter+1
    temp_y=Y0-lam*f.deriv(Y0)
    X1=box_threshold(temp_y)
    t1=(1+sqrt(1+4*t0^2))/2
    Y1=X1+(t0-1)/t1*(X1-X0)
    dis=sqrt(sum((Y1-Y0)^2))

    Y0=Y1
    X0=X1
    t0=t1
  }
  return(X1)
}
