cv.KPR <- function(X, y, K = 5, lambda = seq(0.001,0.2,0.005),
                   kernel = c("matern","exponential")[1],
                   nu = 4.5, power = 1.95, rho = seq(0.05,0.5,0.05), ...){
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if(is.factor(y)){
    y <- as.factor(y)
  }
  
  para.mx <- expand.grid(rho, lambda)[,2:1]
  n <- length(y)
  all.folds <- split(sample(1:n), rep(1:K, length = n))
  pmse.mx <- matrix(0, K, nrow(para.mx))
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    for(j in 1:nrow(para.mx)){
      phat <- try(KPR(X[-omit, , drop = FALSE], y[-omit], X[omit, , drop = FALSE],
                      lambda = para.mx[j,1], kernel = kernel, nu = nu, power = power,
                      rho = para.mx[j,2], ...),  silent = TRUE)
      if(class(phat)=="try-error"){
        pmse.mx[i,j] <- NA
      }else{
        pmse.mx[i,j] <- mean((as.numeric(y[omit])-phat)^2)
      }
    }
  }
  
  MSE.vt <- apply(pmse.mx,2,mean)

  if(all(is.na(MSE.vt))){
    stop("Try other lambda and rho.")
  }else{
    return(list(lambda = para.mx[which.min(MSE.vt),1],
                rho = para.mx[which.min(MSE.vt),2]))
  }
}
