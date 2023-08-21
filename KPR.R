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

#' Kernel Poisson regression
#'
#' Learns a kernel Poisson regression model for count data
#'
#' The method operates by constructing iteratively re-weighted least squares approximations
#' of the log-likelihood loss function and then calling the kernel ridge regression routine
#' to solve those approximations. The least squares approximations are obtained via the Taylor series
#' expansion about the current parameter estimates.
#' 
#' @param K n-by-n matrix of pairwise kernel values over a set of n samples
#' @param y n-by-1 vector of binary response labels
#' @param lambda scalar, regularization parameter
#' @param max.iter maximum number of iterations
#' @param eps convergence precision
#' @param v.init initial parameter estimate for the kernel weights
#' @param b.init initial parameter estimate for the bias term
#' @param silent set to TRUE to suppress run-time output to stdout (default: FALSE)
#' @param balanced boolean specifying whether the balanced model is being evaluated
#' @return A list with two elements:
#' \describe{
#'   \item{v}{n-by-1 vector of kernel weights}
#'   \item{b}{scalar, bias term for the model}
#' }
#' @seealso \code{\link{gelnet.krr}}, \code{\link{gelnet.logreg.obj}}
#' @noRd
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

#' Kernel Poisson Regression
#'
#' @description The function performs a kernel logistic regression for binary outputs.
#'
#' @param X a design matrix with dimension \code{n} by \code{d}, the range of whose values is from 0 to 1.
#' @param y a response vector with length \code{n}. The values in the vector are 0 or 1.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param lambda a positive value specifing the tuning parameter for KLR. The default is 0.01.
#' @param kernel "matern" or "exponential" which specifies the matern kernel or power exponential kernel. The default is "matern".
#' @param nu a positive value specifying the order of matern kernel if \code{kernel} == "matern". The default is 2.5 if matern kernel is chosen.
#' @param power a positive value (between 1.0 and 2.0) specifying the power of power exponential kernel if \code{kernel} == "exponential". The default is 1.95 if power exponential kernel is chosen.
#' @param rho a positive value specifying the scale parameter of matern and power exponential kernels. The default is 0.1.
#'
#' @details This function performs a kernel logistic regression, where the kernel can be assigned to Matern kernel or power exponential kernel by the argument \code{kernel}. The arguments \code{power} and \code{rho} are the tuning parameters in the power exponential kernel function, and \code{nu} and \code{rho} are the tuning parameters in the Matern kernel function. The power exponential kernel has the form \deqn{K_{ij}=\exp(-\frac{\sum_{k}{|x_{ik}-x_{jk}|^{power}}}{rho}),} and the Matern kernel has the form \deqn{K_{ij}=\prod_{k}\frac{1}{\Gamma(nu)2^{nu-1}}(2\sqrt{nu}\frac{|x_{ik}-x_{jk}|}{rho})^{nu} \kappa(2\sqrt{nu}\frac{|x_{ik}-x_{jk}|}{rho}).} The argument \code{lambda} is the tuning parameter for the function smoothness.
#'
#' @return Predictive probabilities at given locations \code{xnew}.
#'
#' @seealso \code{\link{cv.KLR}} for performing cross-validation to choose the tuning parameters.
#' @references Zhu, J. and Hastie, T. (2005). Kernel logistic regression and the import vector machine. Journal of Computational and Graphical Statistics, 14(1), 185-205.
#'
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import GPfit
#' @import gelnet
#'
#' @examples
#' library(calibrateBinary)
#'
#' set.seed(1)
#' np <- 10
#' xp <- seq(0,1,length.out = np)
#' eta_fun <- function(x) exp(exp(-0.5*x)*cos(3.5*pi*x)-1) # true probability function
#' eta_x <- eta_fun(xp)
#' yp <- rep(0,np)
#' for(i in 1:np) yp[i] <- rbinom(1,1, eta_x[i])
#'
#' x.test <- seq(0,1,0.001)
#' etahat <- KLR(xp,yp,x.test)
#'
#' plot(xp,yp)
#' curve(eta_fun, col = "blue", lty = 2, add = TRUE)
#' lines(x.test, etahat, col = 2)
#'
#' #####   cross-validation with K=5    #####
#' ##### to determine the parameter rho #####
#'
#' cv.out <- cv.KLR(xp,yp,K=5)
#' print(cv.out)
#'
#' etahat.cv <- KLR(xp,yp,x.test,lambda=cv.out$lambda,rho=cv.out$rho)
#'
#' plot(xp,yp)
#' curve(eta_fun, col = "blue", lty = 2, add = TRUE)
#' lines(x.test, etahat, col = 2)
#' lines(x.test, etahat.cv, col = 3)
#'
#' @export

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
  #return(sum(diag(solve( K1 %*% A %*% t(K1) + lambda * n * K0, K1))))
  return(sum(diag(t(K1) %*% solve( K1 %*% A %*% t(K1) + lambda * n * K0) %*% K1%*% A)))
}

dispersion <- function(y, X, lambda, kernel, 
                       nu = 4.5, power = 1.95, rho = 0.1){
  
  yhat <- KPR(X,y,X,lambda=lambda,rho=rho, b.init = log(mean(y)), kernel=kernel, nu=nu, power=power)
  D <- 2*sum((ifelse(y==0,0,y*log(y/yhat)) - (y-yhat)))
  D.edf <- edf(X=X, lambda=lambda, kernel=kernel, nu=nu, rho=rho, a=yhat)
  return(list(D=D, D.edf=D.edf, disp=D/D.edf))
}
