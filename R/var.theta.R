var.theta <- function(np, eta, f, x, theta, method=c("L2","LS")[1], emulate.fg=FALSE, ...){
  
  if(is.null(nrow(x))) x <- matrix(x,ncol=1)
  
  n <- nrow(x)
  q <- length(theta)
  
  if(!emulate.fg){  # true simulator (when simulator is cheap)
    df <- jacobian(function(para) {f(x, para)}, theta, ...)
    W <- t(df) %*% diag(eta, n) %*% df / n
    if(method=="LS") W <- W + c(t(df) %*% diag((eta-df)^2, n) %*% df) / n
    V <- matrix(0,q,q)
    for(j in 1:q){
      V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1))^2, para2, ...)[,j], theta, ...)
      V[j,] <- apply(V.tmp, 2, mean)
    }
  }else{ # emulator that approximates the simulator (when simulator is expensive)
    df <- jacobian(function(para) {f(x, para)$mean}, theta, ...)
    W <- t(df) %*% diag(eta, n) %*% df / n
    if(method=="LS") W <- W + c(t(df) %*% diag((eta-df)^2, n) %*% df) / n
    V <- matrix(0,q,q)
    for(j in 1:q){
      V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1)$mean)^2 + f(x, para1)$s2, para2, ...)[,j], theta, ...)
      V[j,] <- apply(V.tmp, 2, mean)
    }
  }
  
  V.inv <- solve(V)
  
  return((4/np)*V.inv%*%W%*%V.inv)
}