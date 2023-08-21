gelnet.krr <- function( K, y, lambda, a, fix.bias=FALSE )
{
  ## Argument verification
  n <- nrow(K)
  stopifnot( length(y) == n )
  stopifnot( length(a) == n )
  stopifnot( is.null(dim(a)) )	## a must not be in matrix form
  
  ## Set up the sample weights and kernalization of the bias term
  A <- diag(a)
  
  if( fix.bias )
  {
    ## Solve the optimization problem in closed form
    m <- solve( K %*% A %*% t(K) + lambda * n * K, K %*% A %*% y )
    list( v = drop(m), b = 0 )
  }
  else
  {
    ## Set up kernalization of the bias term
    K1 <- rbind( K, 1 )
    K0 <- cbind( rbind( K, 0 ), 0 )
    
    ## Solve the optimization problem in closed form
    m <- solve( K1 %*% A %*% t(K1) + lambda * n * K0, K1 %*% A %*% y )
    
    ## Retrieve the weights and the bias term
    list( v = m[1:n], b = m[n+1] )
  }
}

gelnet.kpr <- function( K, y, lambda, max.iter = 100, eps = 1e-5,
                        v.init=rep(0,nrow(K)), b.init=0, silent=FALSE)
{
  ## Argument verification
  stopifnot( nrow(K) == ncol(K) )
  stopifnot( nrow(K) == length(y) )
  
  ## Objective function to minimize
  fobj <- function( v, b )
  {
    s <- K %*% v + b
    R2 <- lambda * t(v) %*% K %*% v / 2
    LL <- mean( y * s - exp(s) )
    R2 - LL
  }
  
  ## Compute the initial objective value
  v <- v.init
  b <- b.init
  fprev <- fobj( v, b )
  
  ## Reduce kernel logistic regression to kernel regression
  ##  about the current Taylor method estimates
  for( iter in 1:max.iter )
  {
    if( silent == FALSE )
    {
      cat( "Iteration", iter, ": " )
      cat( "f =", fprev, "\n" )
    }
    
    ## Compute the current fit
    s <- drop(K %*% v + b)
    
    ## Compute the sample weights and the response
    a <- exp(s)
    z <- s + (y - exp(s)) / a
    
    ## Run the coordinate descent for the resulting regression problem
    m <- gelnet.krr( K, z, lambda, a )
    v <- m$v
    b <- m$b
    
    ## Evaluate the objective function and check convergence criteria
    f <- fobj( v, b )
    if( abs(f - fprev) / abs(fprev) < eps ) break
    else fprev <- f
  }
  
  list( v = v, b = b )
}

# kernel Poisson regression
KPR <- function(X, y, xnew, lambda = 0.01,
                kernel = c("matern","exponential")[1],
                nu = 4.5, power = 1.95, rho = 0.1, ...){
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (any(y < 0)) {
    stop("y needs to be a positive integer.")
  }
  if (!is.matrix(xnew)) {
    xnew <- as.matrix(xnew)
  }
  
  # computer K
  beta <- -log10(rep(rho,ncol(X)))
  if(kernel == "matern") {
    K <- corr_matrix(X, beta, corr = list(type = "matern", nu = nu))
  }else{
    K <- corr_matrix(X, beta, corr = list(type = "exponential", power = power))
  }
  fit <- gelnet.kpr(K, y, lambda, silent = TRUE, ...)
  n <- nrow(X)
  
  phi_new <- apply(xnew, 1, function(xn){
    xn <- matrix(xn, nrow = 1)
    if(kernel == "matern") {
      temp <- 10^beta
      temp <- matrix(temp, ncol = ncol(xn), nrow = (length(X)/ncol(xn)),
                     byrow = TRUE)
      temp <- 2 * sqrt(nu) * abs(X - as.matrix(rep(1, n)) %*%
                                   (xn)) * (temp)
      ID <- which(temp == 0)
      rd <- (1/(gamma(nu) * 2^(nu - 1))) * (temp^nu) * besselK(temp,
                                                               nu)
      rd[ID] <- 1
      return(matrix(apply(rd, 1, prod), ncol = 1))
    }else{
      return(exp(-(abs(X - as.matrix(rep(1, n)) %*% (xn))^power) %*%
                   (10^beta)))
    }
  })
  
  f <- c(fit$v %*% phi_new)+fit$b
  
  return(exp(f))
}


# computing the effective d.o.f.
edf <- function(X, lambda, kernel, 
                nu = 4.5, power = 1.95, rho = 0.1, 
                a){
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  n <- nrow(X)
  
  # computer K
  beta <- -log10(rep(rho,ncol(X)))
  if(kernel == "matern") {
    K <- corr_matrix(X, beta, corr = list(type = "matern", nu = nu))
  }else{
    K <- corr_matrix(X, beta, corr = list(type = "exponential", power = power))
  }
  
  A <- diag(a)
  ## Set up kernalization of the bias term
  K1 <- rbind( K, 1 )
  K0 <- cbind( rbind( K, 0 ), 0 )
  
  ## compute edf
  return(sum(diag(t(K1) %*% solve( K1 %*% A %*% t(K1) + lambda * n * K0) %*% K1%*% A)))
}

# dispersion function
dispersion <- function(y, X, lambda, kernel, 
                       nu = 4.5, power = 1.95, rho = 0.1){
  
  yhat <- KPR(X,y,X,lambda=lambda,rho=rho, b.init = log(mean(y)), kernel=kernel, nu=nu, power=power)
  D <- 2*sum((ifelse(y==0,0,y*log(y/yhat)) - (y-yhat)))
  D.edf <- edf(X=X, lambda=lambda, kernel=kernel, nu=nu, rho=rho, a=yhat)
  return(list(D=D, D.edf=D.edf, disp=D/D.edf))
}
