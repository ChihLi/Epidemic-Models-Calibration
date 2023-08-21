emulation.fit <- function(n, d1, d2, n.rep, lower, upper, 
                          f_xtheta, K = 5, nvar.max = 300, parallel = FALSE, 
                          X = NULL, Y = NULL, 
                          method = c("MRFA", "hetGP")[1], ...){
  
  if(is.null(X) & is.null(Y)){
    set.seed(1)
    X <- maximinLHS(n = n, k = d1 + d2)
    Y <- rep(0, n*n.rep)
    
    X <- t(t(X) * (upper - lower)  + lower)
    
    for(i in 1:n){
      lambda <- f_xtheta(X[i,1:d1], X[i,(d1+1):(d1+d2)])
      Y[((i-1)*n.rep+1) :(i*n.rep)] <- rpois(n.rep, rep(lambda, n.rep))
    }
    X <- matrix(rep(X, each = n.rep), ncol = ncol(X))
  }
  if(method == "MRFA"){
    MRFA.fit <- MRFA_fit(X, Y, converge.tol = 1e-8, nvar.max = nvar.max, model = PoissReg(), parallel = parallel)
    cv.out <- cv.MRFA(X, Y, K = K, lambda = MRFA.fit$lambda, model = PoissReg(), plot.it = FALSE, parallel = parallel)
    lambda_cv <- cv.out$lambda[which.min(cv.out$cv)]
    
    return(list(model.fit = MRFA.fit, lambda_cv = lambda_cv, parallel = parallel, method = method))
  }else{
    hetGP.model <- mleHetGP(X = X, Z = Y, ...)
    return(list(model.fit = hetGP.model, method = method))
  }
}

emulation.pred <- function(x, fit, s2.fit=FALSE)   {
  if(fit$method == "MRFA"){
    PoissReg()@invlink(predict.MRFA(fit$model.fit, x, lambda = fit$lambda_cv, parallel = fit$parallel)$y_hat)
  }else{
    if(s2.fit){
      pred.out <- predict(fit$model.fit, x)
      return(list(mean=pred.out$mean,
                  s2=pred.out$sd2))
    }else{
      predict(fit$model.fit, x)$mean
    }
  }
}