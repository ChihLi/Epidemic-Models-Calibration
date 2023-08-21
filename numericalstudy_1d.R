library(GPfit)
library(gelnet)
library(ggplot2)
library(plyr)
library(lhs)
library(scales)
library(MRFA)
library(grplasso)
library(foreach)
library(doParallel)
library(snowfall)
library(hetGP)
library(numDeriv)
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
np <- 50 # sample size of physical data
xp <- seq(0,2*pi,length.out = np)
eta_fun <- function(x) exp(x/2)*sin(x/2) + 30 # true process
eta_x <- eta_fun(xp)

# simulation model
f_xtheta <- function(x,theta) {
  eta_fun(x) - 5*sqrt(theta^2-theta+1) * (sin(theta * x) + cos(theta * x))
}

#true L2 theta
true.cpara <- -0.1789

# emualtion comparison 
sample.size <- matrix(c(25, 25, 50, 100, 50, 100, 50, 100), ncol = 2)
training.time <- test.time <- mse <- matrix(0, nrow = nrow(sample.size), ncol = 2)
for(i in 1:nrow(sample.size)){
  # MRFA
  start.time <- Sys.time()
  MRFA.fit <- emulation.fit(n = sample.size[i,1], d1 = 1 , d2 = 1, n.rep = sample.size[i,2], 
                            nvar.max = 300, K = 10, parallel = TRUE,
                            lower = c(0,-1), upper = c(2*pi, 1), 
                            f_xtheta = f_xtheta, method = "MRFA")
  end.time <- Sys.time()
  training.time[i,1] <- difftime(end.time, start.time, units = "secs")
  
  # hetGP
  start.time <- Sys.time()
  hetGP.fit <- emulation.fit(n = sample.size[i,1], d1 = 1 , d2 = 1, n.rep = sample.size[i,2], 
                             lower = c(0,-1), upper = c(2*pi, 1),  f_xtheta = f_xtheta, 
                             method = "hetGP")
  end.time <- Sys.time()
  training.time[i,2] <- difftime(end.time, start.time, units = "secs")
  
  # test the emulation performance
  set.seed(i)
  xtest <- matrix(runif(10000*2), ncol = 2)
  xtest[,1] <- xtest[,1] * 2*pi
  xtest[,2] <- (xtest[,2] - 0.5)*2
  start.time <- Sys.time()
  mse[i,1] <- mean((emulation.pred(xtest, MRFA.fit) - f_xtheta(xtest[,1], xtest[,2]))^2)
  end.time <- Sys.time()
  test.time[i,1] <- difftime(end.time, start.time, units = "secs")
  start.time <- Sys.time()
  mse[i,2] <- mean((emulation.pred(xtest, hetGP.fit) - f_xtheta(xtest[,1], xtest[,2]))^2)
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
lower <- -1
upper <- 1
hetGP.fit <- emulation.fit(n = 25, d1 = 1 , d2 = 1, n.rep = 50, 
                           lower = c(0,lower), upper = c(2*pi, upper),  f_xtheta = f_xtheta, 
                           method = "hetGP")

# test input
x.test <- seq(0,2*pi,length.out = 100)

# run calibration for 100 repetitions
l2.par <- l2.emu.par <- ols.par <- ols.emu.par <- mle.par <- mle.emu.par <- 
  l2.var <- l2.emu.var <- rep(0,100)
for(ii in 1:100){
  set.seed(ii)
  # generate the physical output
  yp <- rep(0,np)
  for(i in 1:np) yp[i] <- rpois(1,eta_x[i])
  
  ## kernel Poisson regression to get the estimated mean process
  cv.out <- cv.KPR(xp,yp,K=10, lambda = seq(0.05, 0.5, 0.05), 
                   rho = 2*pi*seq(0.5, 3, 0.5), b.init = log(mean(yp)))
  etahat.cv <- KPR(xp,yp,x.test,lambda=cv.out$lambda,rho=cv.out$rho, b.init = log(mean(yp)))
  
  ## L2 calibration: without emulation
  ini.val <- true.cpara
  L2_fun <- function(theta) {
    mean((etahat.cv - f_xtheta(x.test, theta))^2)
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], L2_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  l2.par[ii] <- opt.sol[which.min(opt.val)]
  l2.var[ii] <- drop(var.theta(np, etahat.cv, f_xtheta, x.test, l2.par[ii], method="L2", emulate.fg=FALSE))
  
  ## L2 calibration: with emulator
  L2_emulation_fun <- function(theta) {
    xnew <- cbind(x.test, matrix(rep(theta, each = length(x.test)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit, s2.fit = TRUE)
    mean((etahat.cv - predictions$mean)^2)+mean(predictions$s2)
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], L2_emulation_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  l2.emu.par[ii] <- opt.sol[which.min(opt.val)]  
  f.emulator <- function(x, theta){
    xnew <- cbind(x, matrix(rep(theta, each = length(x)), ncol = length(theta)))
    return(emulation.pred(x = xnew, fit = hetGP.fit, s2.fit = TRUE))
  }
  l2.emu.var[ii] <- drop(var.theta(np, etahat.cv, f.emulator, x.test, l2.emu.par[ii], method="L2", emulate.fg=TRUE))
  
  ## Least squares approach: without emulation
  OLS_fun <- function(theta) {
    mean((yp - f_xtheta(xp, theta))^2)
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], OLS_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  ols.par[ii] <- opt.sol[which.min(opt.val)]

  ## Least squares approach: with emulation
  OLS_emulation_fun <- function(theta) {
    xnew <- cbind(xp, matrix(rep(theta, each = length(xp)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit, s2.fit = TRUE)
    mean((yp - predictions$mean)^2) + mean(predictions$s2)
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], OLS_emulation_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  ols.emu.par[ii] <- opt.sol[which.min(opt.val)]

  ## MLE: without emulation
  MLE_fun <- function(theta) {
    -mean((yp*log(f_xtheta(xp, theta)) - f_xtheta(xp, theta)))
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], MLE_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  mle.par[ii] <- opt.sol[which.min(opt.val)]
  
  ## MLE: with emulation
  MLE_emulation_fun <- function(theta) {
    xnew <- cbind(xp, matrix(rep(theta, each = length(xp)), ncol = length(theta)))
    predictions <- emulation.pred(x = xnew, fit =  hetGP.fit)
    predictions[predictions<=0] <- 1e-8
    -mean((yp*log(predictions) - predictions))
  }
  opt.out <- vector("list", length(ini.val))
  for(i in 1:length(ini.val)) opt.out[[i]] <- optim(ini.val[i], MLE_emulation_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- sapply(opt.out, function(x) x$par)
  mle.emu.par[ii] <- opt.sol[which.min(opt.val)]
  
}

##  build confidence interval
ub <- l2.par + qnorm(0.975) * sqrt(l2.var)
lb <- l2.par - qnorm(0.975) * sqrt(l2.var)
cat("coverage of 0.95 confidence interval without emulator:", mean(ub > true.cpara & lb < true.cpara), "\n")

ub <- l2.emu.par + qnorm(0.975) * sqrt(l2.emu.var)
lb <- l2.emu.par - qnorm(0.975) * sqrt(l2.emu.var)
cat("coverage of 0.95 confidence interval with emulator:", mean(ub > true.cpara & lb < true.cpara), "\n")

## visualize the MSE
MSE <- c(mean((l2.par - true.cpara)^2), 
         mean((ols.par - true.cpara)^2), 
         mean((mle.par - true.cpara)^2),
         mean((l2.emu.par - true.cpara)^2), 
         mean((ols.emu.par - true.cpara)^2), 
         mean((mle.emu.par - true.cpara)^2))

df <- data.frame(emulator=c(rep("true",3), rep("emulator",3)),
                 method=c("L2", "OLS", "MLE", "L2+emulator", "OLS+emulator", "MLE+emulator"),
                 MSE=MSE)
df[,"method"]  <- factor(df[,"method"] , levels = c("L2", "OLS", "MLE", "L2+emulator", "OLS+emulator", "MLE+emulator"))

df[,"q1"] <- c(quantile((l2.par - true.cpara)^2, 0.05), 
                quantile((ols.par - true.cpara)^2, 0.05), 
                quantile((mle.par - true.cpara)^2, 0.05),
                quantile((l2.emu.par - true.cpara)^2, 0.05), 
                quantile((ols.emu.par - true.cpara)^2, 0.05), 
                quantile((mle.emu.par - true.cpara)^2, 0.05))
df[,"q3"] <- c(quantile((l2.par - true.cpara)^2, 0.95), 
                quantile((ols.par - true.cpara)^2, 0.95), 
                quantile((mle.par - true.cpara)^2, 0.95),
                quantile((l2.emu.par - true.cpara)^2, 0.95), 
                quantile((ols.emu.par - true.cpara)^2, 0.95), 
                quantile((mle.emu.par - true.cpara)^2, 0.95))

p.mse <-ggplot(df, aes(x=method, y=MSE, fill=emulator)) + ylab("MSE") + xlab("") +
  geom_bar(stat="identity")+theme_minimal() + theme(legend.position="none") + 
  geom_text(aes(y=MSE, label=round(MSE,5)), vjust=-0.3, color="black", size=3.5) + 
  scale_y_continuous(limits = c(0, max(df$q3)), oob = rescale_none)+
  scale_x_discrete(labels=c("L2" = "L2", "OLS" = "LSE", "MLE" = "MLE", "L2+emulator" = "L2+\nemulator", "OLS+emulator" = "LSE+\nemulator", "MLE+emulator" = "MLE+\nemulator")) 

# Figure 1 (right)
p.mse <- p.mse + geom_errorbar(aes(x=method, ymin=q1, ymax=q3), width=0.2, alpha=0.5)

# Figure 1 (left)
set.seed(100)
yp <- rep(0,np)
for(i in 1:np) yp[i] <- rpois(1,eta_x[i])
p.demo <- ggplot(data.frame(x=xp, y=yp), aes(x,y)) + geom_point() +
  stat_function(fun=function(x) eta_fun(x), aes(linetype = "solid") ) +
  stat_function(fun=function(x) f_xtheta(x, true.cpara), aes(linetype = "dashed") ) +
  scale_linetype_manual("", values = c("dashed","solid"),
                        labels = c(expression(italic(f(x,theta^"*"))),expression(italic(lambda(x))))) +
  theme_minimal() + theme(legend.position = "top", legend.text = element_text(size=14)) + labs(fill = "")

## Figure 1
# save the plot 
pdf("example_1d.pdf", width = 9, height = 3.5)
multiplot(p.demo,p.mse, layout = matrix(c(1,2),nrow=1,ncol=2))
dev.off()

