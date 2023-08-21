library(GPfit)
library(gelnet)
library(ggplot2)
library(plyr)
library(lhs)
library(MRFA)
library(scales)
library(grplasso)
library(foreach)
library(doParallel)
library(snowfall)
library(hetGP)
library(numDeriv)
library(plyr)
source("R/emulation.fun.R")
source("R/var.theta.R")
source("R/multiplot.R")
source("R/cv.KPR.R")
source("R/KPR.R")

# parallel computing for emulator comparison
sfInit(parallel = TRUE, cpus = 10L)
cl <- sfGetCluster()
registerDoParallel(cl)

set.seed(1)
np <- 50  # sample size of physical data
xp <- seq(0,2,length.out = np)
eta_fun <- function(x) 3*x + 3*x*sin(5*x) +3 # true process
eta_x <- eta_fun(xp)

# simulation model
f_xtheta <- function(x,theta) {
  theta[1] + theta[2]*x + theta[3]*x^2
}

#true L2 theta
true.cpara <- optim(c(0,0,0), fn = function(g) {
  x.grid <- seq(0,2,length.out = 1000)
  mean((eta_fun(x.grid) - f_xtheta(x.grid, g))^2)
},
lower = 0, upper = 5, method = "L-BFGS-B")$par

# emualtion comparison 
sample.size <- matrix(c(300, 300, 500, 500, 500, 50, 100, 5, 25, 50), ncol = 2)
training.time <- test.time <- mse <- matrix(0, nrow = nrow(sample.size), ncol = 2)
for(i in 1:nrow(sample.size)){
  # MRFA
  start.time <- Sys.time()
  MRFA.fit <- emulation.fit(n = sample.size[i,1], d1 = 1 , d2 = 3, n.rep = sample.size[i,2], 
                            nvar.max = 300, K = 10, parallel = TRUE,
                            lower = c(0,0,0,0), upper = c(2, 5, 5, 5), 
                            f_xtheta = f_xtheta, method = "MRFA")
  end.time <- Sys.time()
  training.time[i,1] <- difftime(end.time, start.time, units = "secs")
  
  # hetGP
  start.time <- Sys.time()
  hetGP.fit <- emulation.fit(n = sample.size[i,1], d1 = 1 , d2 = 3, n.rep = sample.size[i,2], 
                             lower = c(0,0,0,0), upper = c(2, 5, 5, 5),  f_xtheta = f_xtheta, 
                             method = "hetGP")
  end.time <- Sys.time()
  training.time[i,2] <- difftime(end.time, start.time, units = "secs")
  
  # test the emulation performance
  set.seed(i)
  xtest <- matrix(runif(10000*4), ncol = 4)
  xtest[,1] <- xtest[,1] * 2
  xtest[,2:4] <- xtest[,2:4] * 5
  actual.sim <- apply(xtest, 1, function(x) f_xtheta(x[1], x[2:4]))
  start.time <- Sys.time()
  mse[i,1] <- mean((emulation.pred(xtest, MRFA.fit) - actual.sim)^2)
  end.time <- Sys.time()
  test.time[i,1] <- difftime(end.time, start.time, units = "secs")
  start.time <- Sys.time()
  mse[i,2] <- mean((emulation.pred(xtest, hetGP.fit) - actual.sim)^2)
  end.time <- Sys.time()
  test.time[i,2] <- difftime(end.time, start.time, units = "secs")
}
sfStop()
rm(cl)

cat("RMSE:\n")
print(sqrt(mse))
cat("Training time:\n")
print(training.time)
cat("Prediction time:\n")
print(test.time)

# training emulator
lower <- 0
upper <- 5
hetGP.fit <- emulation.fit(n = 300, d1 = 1 , d2 = 3, n.rep = 100, 
                           lower = c(0,lower,lower,lower), upper = c(2, upper, upper, upper),  f_xtheta = f_xtheta, 
                           method = "hetGP")

# test input
x.test <- seq(0,2,length.out = 100)

# run calibration for 100 repetitions
l2.par <- l2.emu.par <- l2.var <- ols.par <- ols.emu.par <- mle.par <- mle.emu.par <- matrix(0,nrow=100,ncol=3)
for(ii in 1:100){
  
  # generate the physical output
  yp <- rep(0,np)
  for(i in 1:np) yp[i] <- rpois(1,eta_x[i])
  
  ## kernel Poisson regression to get the estimated mean process
  cv.out <- cv.KPR(xp,yp,K=10, lambda = seq(0.05, 0.5, 0.05), 
                   rho = 2*seq(0.5, 5, 0.5), b.init = log(mean(yp)))
  etahat.cv <- KPR(xp,yp,x.test,lambda=cv.out$lambda,rho=cv.out$rho, b.init =log(mean(yp)))
  
  ## L2 calibration: without emulation
  L2_fun <- function(theta) {
    mean((etahat.cv - f_xtheta(x.test, theta))^2)
  }
  ini.mx <- expand.grid(c(1.5,3.5),c(1.5,3.5),c(1.5,3.5))
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], L2_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  l2.par[ii,] <- opt.sol[[which.min(opt.val)]]
  l2.var[ii,] <- diag(var.theta(np, etahat.cv, f_xtheta, x.test, l2.par[ii,], method="L2"))

  ## L2 calibration: with emulator
  L2_emulation_fun <- function(theta) {
    xnew <- cbind(x.test, matrix(rep(theta, each = length(x.test)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit, s2.fit = TRUE)
    mean((etahat.cv - predictions$mean)^2)+ mean(predictions$s2)
  }
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], L2_emulation_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  l2.emu.par[ii,] <- opt.sol[[which.min(opt.val)]]

  ## Least squares approach: without emulation
  OLS_fun <- function(theta) {
    mean((yp - f_xtheta(xp, theta))^2)
  }
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], OLS_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  ols.par[ii,] <- opt.sol[[which.min(opt.val)]]
  
  ## Least squares approach: with emulation
  OLS_emulation_fun <- function(theta) {
    xnew <- cbind(xp, matrix(rep(theta, each = length(xp)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit, s2.fit = TRUE)
    mean((yp - predictions$mean)^2) + mean(predictions$s2)
  }
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], OLS_emulation_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  ols.emu.par[ii,] <- opt.sol[[which.min(opt.val)]]
  
  ## MLE: without emulation
  MLE_fun <- function(theta) {
    -mean((yp*log(f_xtheta(xp, theta)) - f_xtheta(xp, theta)))
  }
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], MLE_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  mle.par[ii,] <- opt.sol[[which.min(opt.val)]]
  
  ## MLE: with emulation
  MLE_emulation_fun <- function(theta) {
    xnew <- cbind(xp, matrix(rep(theta, each = length(xp)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit)
    predictions[predictions<=0] <- 1e-8
    -mean((yp*log(predictions) - predictions))
  }
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], MLE_emulation_fun, lower = rep(lower,3), upper = rep(upper,3), method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  mle.emu.par[ii,] <- opt.sol[[which.min(opt.val)]]
}

## visualize the MSE
mse <- function(x) (x-true.cpara)^2
MSE <- rbind(apply(apply(l2.par, 1, mse),1,mean), 
             apply(apply(ols.par, 1, mse),1,mean),
             apply(apply(mle.par, 1, mse),1,mean),
             apply(apply(l2.emu.par, 1, mse),1,mean),
             apply(apply(ols.emu.par, 1, mse),1,mean),
             apply(apply(mle.emu.par, 1, mse),1,mean))

df <- data.frame(emulator=c(rep("true",3), rep("emulator",3)),
                 method=c("L2", "OLS", "MLE", "L2+emulator", "OLS+emulator", "MLE+emulator"),
                 MSE=MSE)
df[,"method"]  <- factor(df[,"method"] , levels = c("L2", "OLS", "MLE", "L2+emulator", "OLS+emulator", "MLE+emulator"))

q1 <- rbind(apply(apply(l2.par, 1, mse),1,quantile,0.05), 
                   apply(apply(ols.par, 1, mse),1,quantile,0.05),
                   apply(apply(mle.par, 1, mse),1,quantile,0.05),
                   apply(apply(l2.emu.par, 1, mse),1,quantile,0.05),
                   apply(apply(ols.emu.par, 1, mse),1,quantile,0.05),
                   apply(apply(mle.emu.par, 1, mse),1,quantile,0.05))
colnames(q1) <- paste0("q1.",1:3)
q3 <- rbind(apply(apply(l2.par, 1, mse),1,quantile,0.95), 
                   apply(apply(ols.par, 1, mse),1,quantile,0.95),
                   apply(apply(mle.par, 1, mse),1,quantile,0.95),
                   apply(apply(l2.emu.par, 1, mse),1,quantile,0.95),
                   apply(apply(ols.emu.par, 1, mse),1,quantile,0.95),
                   apply(apply(mle.emu.par, 1, mse),1,quantile,0.95))
colnames(q3) <- paste0("q3.",1:3)
df <- cbind(df, q1, q3)


p1.mse <-ggplot(df, aes(x=method, y=MSE.1, fill=emulator)) + ylab("MSE") + xlab(expression(theta[1])) +
  geom_bar(stat="identity")+theme_minimal() + theme(legend.position="none") + 
  geom_text(aes(y=MSE.1, label=round(MSE.1,4)), vjust=-0.3, color="black", size=3.5) + 
  scale_y_continuous(limits = c(min(df$q1.1), max(df$q3.1)), oob = rescale_none)+
  scale_x_discrete(labels=c("L2" = "L2", "OLS" = "LSE", "MLE" = "MLE", "L2+emulator" = "L2+\nemulator", "OLS+emulator" = "LSE+\nemulator", "MLE+emulator" = "MLE+\nemulator")) 
p1.mse <- p1.mse + geom_errorbar(aes(x=method, ymin=q1.1, ymax=q3.1), width=0.2, alpha=0.5)


p2.mse <-ggplot(df, aes(x=method, y=MSE.2, fill=emulator)) + ylab("MSE") + xlab(expression(theta[2])) +
  geom_bar(stat="identity")+theme_minimal() + theme(legend.position="none") + 
  geom_text(aes(y=MSE.2, label=round(MSE.2,4)), vjust=-0.3, color="black", size=3.5) + 
  scale_y_continuous(limits = c(min(df$q1.2), max(df$q3.2)), oob = rescale_none)+
  scale_x_discrete(labels=c("L2" = "L2", "OLS" = "LSE", "MLE" = "MLE", "L2+emulator" = "L2+\nemulator", "OLS+emulator" = "LSE+\nemulator", "MLE+emulator" = "MLE+\nemulator")) 
p2.mse <- p2.mse + geom_errorbar(aes(x=method, ymin=q1.2, ymax=q3.2), width=0.2, alpha=0.5)


p3.mse <-ggplot(df, aes(x=method, y=MSE.3, fill=emulator)) + ylab("MSE") + xlab(expression(theta[3])) +
  geom_bar(stat="identity")+theme_minimal() + theme(legend.position="none") + 
  geom_text(aes(y=MSE.3, label=round(MSE.3,4)), vjust=-0.3, color="black", size=3.5) + 
  scale_y_continuous(limits = c(min(df$q1.3), max(df$q3.3)), oob = rescale_none)+
  scale_x_discrete(labels=c("L2" = "L2", "OLS" = "LSE", "MLE" = "MLE", "L2+emulator" = "L2+\nemulator", "OLS+emulator" = "LSE+\nemulator", "MLE+emulator" = "MLE+\nemulator")) 
p3.mse <- p3.mse + geom_errorbar(aes(x=method, ymin=q1.3, ymax=q3.3), width=0.2, alpha=0.5)

## Figure 2
# save the plot 
pdf("example_3d.pdf", width = 12, height = 3.5)
multiplot(p1.mse, p2.mse, p3.mse, layout = matrix(1:3,ncol=3))
dev.off()

