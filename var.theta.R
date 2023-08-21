var.theta <- function(np, eta, f, x, theta, method=c("L2","LS")[1], emulate.fg=FALSE, ...){
  
  if(is.null(nrow(x))) x <- matrix(x,ncol=1)
  
  n <- nrow(x)
  q <- length(theta)
  
  # if(!emulate.fg){
  #   # true simulator (when simulator is cheap)
  #   p.grad <- function(para, i){
  #     p.tmp <- function(para, a) c(f(a, para))
  #     rootSolve::gradient(p.tmp, para, a = x[i,], centered = FALSE, pert = 1e-4)
  #   }
  # 
  #   L.hessian <- function(para, i){
  #     l.tmp <- function(para, a) c((eta[i] - f(a, para))^2)
  #     rootSolve::hessian(l.tmp, para, a = x[i,], centered = FALSE, pert = 1e-4)
  #   }
  # }else{
  #   # emulator that approximates the simulator (when simulator is expensive)
  #   p.emulator <- function(x, theta){
  #     xnew <- cbind(x, matrix(rep(theta, each = length(x)), ncol = length(theta)))
  #     return(emulation.pred(x = xnew, fit = p, s2.fit = TRUE))
  #   }
  #   p.grad <- function(para, i){
  #     p.tmp <- function(para, a) c(p.emulator(a, para)$mean)
  #     gradient(p.tmp, para, a = x[i,], centered = FALSE, pert = 1e-4)
  #   }
  #   L.hessian <- function(para, i){
  #     l.tmp <- function(para, a) c((eta[i] - p.emulator(a, para)$mean)^2+p.emulator(a, para)$s2)
  #     rootSolve::hessian(l.tmp, para, a = x[i,], centered = FALSE, pert = 1e-4)
  #   }
  # }
  # 
  # W <- W0 <- V <- matrix(0,length(theta),length(theta))
  # for(i in 1:n){
  #   W <- W + (1/n)*eta[i]*t(p.grad(theta, i))%*%p.grad(theta, i)
  #   V <- V + (1/n)*L.hessian(theta, i)
  #   if(method=="LS") W0 <- W0 + (1/n)*(eta[i]-p.grad(theta, i))^2*t(p.grad(theta, i))%*%p.grad(theta, i)
  # }
  # W <- W + W0
  # V.inv <- solve(V)
  
  # if(!emulate.fg){  # true simulator (when simulator is cheap)
  #   df <- jacobian(function(para) {f(x, para)}, theta)
  #   W <- t(df) %*% diag(eta, n) %*% df / n
  #   if(method=="LS") W <- W + c(t(df) %*% diag((eta-df)^2, n) %*% df) / n
  #   V <- matrix(0,q,q)
  #   for(j in 1:q){
  #     V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1))^2, para2)[,j], theta)
  #     V[j,] <- apply(V.tmp, 2, mean)
  #   }
  # }else{ # emulator that approximates the simulator (when simulator is expensive)
  #   df <- jacobian(function(para) {f(x, para)$mean}, theta)
  #   W <- t(df) %*% diag(eta, n) %*% df / n
  #   if(method=="LS") W <- W + c(t(df) %*% diag((eta-df)^2, n) %*% df) / n
  #   V <- matrix(0,q,q)
  #   for(j in 1:q){
  #     V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1)$mean)^2 + f(x, para1)$s2, para2)[,j], theta)
  #     V[j,] <- apply(V.tmp, 2, mean)
  #   }
  # }
  
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
  
  # if(!emulate.fg){  # true simulator (when simulator is cheap)
  # 
  #   df <- jacobian(function(para) {f(x, para)}, theta)
  #   W <- V <- matrix(0,q,q)
  #   for(i in 1:n){
  #     W <- W + eta[i] * t(df[i,,drop=FALSE]) %*% df[i,,drop=FALSE] / n
  #     V <- V + numDeriv::hessian(function(para) (eta[i] - f(x, para)[i])^2, theta)/n
  #   }
  # 
  #   # V <- matrix(0,q,q)
  #   # for(j in 1:q){
  #   #   V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1))^2, para2, method="simple")[,j], theta, method="simple")
  #   #   V[j,] <- apply(V.tmp, 2, mean)
  #   # }
  # }else{ # emulator that approximates the simulator (when simulator is expensive)
  #   df <- jacobian(function(para) {f(x, para)$mean}, theta)
  #   W <- t(df) %*% diag(eta, n) %*% df / n
  #   if(method=="LS") W <- W + c(t(df) %*% diag((eta-df)^2, n) %*% df) / n
  #   V <- matrix(0,q,q)
  #   for(j in 1:q){
  #     V.tmp <- jacobian(function(para2) jacobian(function(para1) (eta - f(x, para1)$mean)^2 + f(x, para1)$s2, para2)[,j], theta)
  #     V[j,] <- apply(V.tmp, 2, mean)
  #   }
  # }
  V.inv <- solve(V)
  
  return((4/np)*V.inv%*%W%*%V.inv)
}