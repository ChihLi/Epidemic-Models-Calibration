library(hetGP)
library(numDeriv)
library(lhs)
library(SimInf)
library(deSolve)
library(GPfit)
library(covid19.analytics)
library(ggplot2)
source("var.theta.R")
source("cv.KPR.R")
source("KPR.R")
source("SEIR.ode.R")
source("SEIR.sim.R")
source("erfinv.R")

### load the reported data
n <- 20 # number of countries we study here 
select.date <- 44:409
all.ts <- covid19.JHU.data(case = "ts-ALL")
all.countries <- as.character(unique(all.ts[,"Country.Region"]))
summary.confirm <- aggregate(all.ts[all.ts[,"status"] == "confirmed", max(select.date)], 
                             list(all.ts[all.ts[,"status"] == "confirmed","Country.Region"]), sum)
country.rank <- summary.confirm[sort.int(summary.confirm$x, decreasing = TRUE, index.return = TRUE)$ix, "Group.1"]

### population of each country
N.world <- read.csv("world_population.csv", stringsAsFactors = FALSE)
N.world[3,1] <- "US"

### select the top n countries
country.select <- country.rank[1:n]

### setup
# estimated parameters (original scale)
para.d <- matrix(0, ncol = 6, nrow = n)
colnames(para.d) <- c("I.init", "E.init", "R.init", "beta", "kappa", "gamma")
para.d <- data.frame("country" = country.select, para.d)
para.s <- para.d
# predictions
pred.d <- uq.d <- lq.d <- pred.s <- uq.s <- lq.s <- 
  data.frame(matrix(0, nrow = length(select.date)-1, ncol = n))
# incubation parameter
incubation <- incubation.sd <- rep(0, 20)
# R0
R0.d <- matrix(0, ncol = 2, nrow = n)
colnames(R0.d) <- c("R0", "sd")
R0.d <- data.frame("country" = country.select, R0.d)
R0.s <- R0.d
# estimated parameters (log scale)
l2.par <- l2.emu.par <- matrix(0,nrow=n,ncol=6)
# overdispersion parameter
phi.hat <- rep(0,n)
# p-value of deviance goodness-of-fit
pval.chi <- rep(0,n)
# daily cases
Infected.out <- data.frame(matrix(0, nrow = length(select.date)-1, ncol = n))
colnames(Infected.out) <- country.select

### running L2 calibration for all the countries
set.seed(1)
for(ii in 1:n){
  country <- country.select[ii]
  cat(country, ":   ")
  
  N <- N.world[N.world[,1] == country,2] # population of the country
  
  # preprocessing
  subdata.ts <- all.ts[all.ts[,"Country.Region"] == country,c(select.date, ncol(all.ts))]
  if(nrow(subdata.ts) > 3){
    subdata.ts <- aggregate(subdata.ts[,-ncol(subdata.ts)], list(subdata.ts[,"status"]), sum)
    subdata.ts <- subdata.ts[,-1]
  }else{
    subdata.ts <- subdata.ts[,-ncol(subdata.ts)]
  }
  subdata.ts <- as.matrix(subdata.ts)
  
  # I and R of the first day
  I.day0 <- subdata.ts[1,1] 
  I.day0 <- max(I.day0, 1)
  R.day0 <- subdata.ts[2,1] + subdata.ts[3,1]
  R.day0 <- max(R.day0, 1)
  daily.case <- c(as.numeric(pmax(0, diff(subdata.ts[1,]))))
  names(daily.case) <- NULL
  Infected.out[,country.select == country] <- daily.case
  yp <- daily.case
  xp <- 1:length(yp)
  xp <- xp[yp>0] # clean misreported zeros
  yp <- yp[yp>0] # clean misreported zeros
  
  # running poisson kernel regression
  x.test <- 1:length(daily.case)
  cv.out <- cv.KPR(xp, yp, K = 5, nu=2.5,
                   lambda = sqrt(median(yp))*seq(0.05, 0.5, 0.05),
                   rho = max(xp)/10*seq(0.5, 3, 0.5), b.init = log(mean(yp)))
  etahat.cv <- KPR(xp,yp,x.test, nu=2.5,lambda=cv.out$lambda,rho=cv.out$rho, b.init = log(mean(yp)))

  # perform GOF and estimate the dispersion parameter
  gof <- dispersion(yp, xp, lambda = cv.out$lambda, 
                            kernel="matern", nu=2.5, 
                            rho = cv.out$rho)
  phi.hat[ii] <- gof$disp
  pval.chi[ii] <- pchisq(gof$D, df=gof$D.edf, lower.tail=FALSE)
  
  ##### deterministic SEIR model #####
  # SEIR model simulation
  f_xtheta <- function(x, theta){
    para <- theta[1:3]
    I.init <- exp(theta[4])
    E.init <- exp(theta[5])
    R.init <- exp(theta[6])
    init <- c(S = N - I.init - E.init - R.init, E = E.init, I = I.init, R = R.init)
    names(para) <- c("beta", "kappa", "gamma")
    out <- ode(y = init, times = x, func = SEIR.ode, parms = para)
    infect.out <- para["kappa"] * out[,"E"]
    infect.out <- infect.out * sum(yp)/sum(infect.out)
    return(infect.out)
  }
  
  # L2 estimates: deterministic SEIR
  L2_fun <- function(theta) {
    sqrt(mean((etahat.cv - f_xtheta(x.test, theta))^2))
  }
  
  # range is from:
  # https://towardsdatascience.com/forecasting-the-covid-19-trend-using-the-seir-model-90979abb9e64
  lower <- c(c(0.05,1/5.6,1/9),log(I.day0*1),log(I.day0*2),log(R.day0*1))
  upper <- c(c(1,1/4.8,1/4),log(I.day0*10),log(I.day0*20),log(R.day0*10))
  
  set.seed(1)
  ini.mx <- maximinLHS(20, 6)
  ini.mx <- t(t(ini.mx) * (upper-lower) + lower)
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], L2_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  l2.par[ii,] <- opt.sol[[which.min(opt.val)]]
  SEIR_par <- l2.par[ii,1:3]
  I.init <- ceiling(exp(l2.par[ii,4]))
  E.init <- ceiling(exp(l2.par[ii,5]))
  R.init <- ceiling(exp(l2.par[ii,6]))
  
  names(SEIR_par) <- c("beta", "kappa", "gamma")
  cat("deterministic: R0 =", SEIR_par[1]/SEIR_par[3], "\n")
  
  para.d[para.d$country == country,-1] <- c(I.init, E.init, R.init, SEIR_par)
  
  # compute variance of R0
  # covariance matrix of parameter estimators
  cov.mx <- var.theta(length(xp), etahat.cv, f_xtheta, x.test, l2.par[ii,], method="L2", side=rep(1,6))
  
  # delta method to compute standard deviation of R0
  R0.d[R0.d[,1] == country, 2] <- SEIR_par[1]/SEIR_par[3]
  R0.d[R0.d[,1] == country, 3] <- sqrt(matrix(c(1/SEIR_par[3],-SEIR_par[1]/(SEIR_par[3])^2), nrow = 1) %*% cov.mx[-c(2,4:6),-c(2,4:6)] %*% matrix(c(1/SEIR_par[3],-SEIR_par[1]/(SEIR_par[3])^2), ncol = 1))
  R0.d[R0.d[,1] == country, 3] <- R0.d[R0.d[,1] == country, 3] * sqrt(phi.hat[ii]) # times dispersion parameter
    
  cat("sd R0 =", R0.d[R0.d[,1] == country, 3], "\n")
  
  # delta method to compute the variance of model outputs
  df <- jacobian(function(para) {f_xtheta(x.test, para)}, l2.par[ii,])
  var.d <- diag(df %*% cov.mx %*% t(df)) * phi.hat[ii]
  # confidence intervals for the model output
  uq.d[,ii] <- f_xtheta(x.test, l2.par[ii,]) + qnorm(0.975)* sqrt(var.d)
  lq.d[,ii] <- f_xtheta(x.test, l2.par[ii,]) - qnorm(0.975)* sqrt(var.d)
  pred.d[,ii] <- f_xtheta(x.test, l2.par[ii,])
  
  ##### stochastic SEIR model #####
  
  # simulate stochastic SEIR model
  # [WARNING] the simulations take a very long time
  # Suggest running these simulations offline
  sim.dat <- SEIR.sim(n.rep=50, n.train=60, N=N,
                      xp=c(xp), 
                      lower=lower, upper=upper)
  
  # fit hetGP model: training emulator
  # (take log for numerical stability)
  hetGP.fit <- mleHetGP(X = sim.dat$Xs, Z=log(sim.dat$Ys+1),
                        covtype = "Matern5_2", settings = list(linkThetas = "none"))
  
  f.emulator <- function(x, theta){
    xnew <- cbind(x, matrix(rep(theta, each = length(x)), ncol = length(theta)))
    predictions <- predict(hetGP.fit, xnew)
    f.mean <- predictions$mean
    f.var <- predictions$sd2
    
    # mean and variance of log normal
    pred.mean <- exp(f.mean-1 + f.var/2)
    pred.var <- (exp(f.var)-1)*exp(2*(f.mean-1) + f.var)
    scale.pred <- sum(yp)/sum(pred.mean)
    pred.mean <- pred.mean*scale.pred
    pred.var <- pred.var*scale.pred^2
    return(list(mean=pred.mean, s2=pred.var, scale.pred=scale.pred))
  }
  
  # L2 estimates: stochastic SEIR
  L2_emulation_fun <- function(theta) {
    predictions <- f.emulator(x.test, theta)
    out <- sqrt(mean((etahat.cv - predictions$mean)^2)+mean(predictions$s2))
    if(is.na(out)) out <- 1e12
    if(!is.finite(out)) out <- 1e12
    return(out)
  }
  
  # [WARNING] this can take a long time
  opt.out <- vector("list", nrow(ini.mx))
  for(i in 1:nrow(ini.mx)) opt.out[[i]] <- optim(ini.mx[i,], L2_emulation_fun, lower = lower, upper = upper, method = "L-BFGS-B")
  opt.val <- sapply(opt.out, function(x) x$value)
  opt.sol <- lapply(opt.out, function(x) x$par)
  l2.emu.par[ii,] <- opt.sol[[which.min(opt.val)]]
  SEIR_par <- l2.emu.par[ii,1:3]
  I.init <- ceiling(exp(l2.emu.par[ii,4]))
  E.init <- ceiling(exp(l2.emu.par[ii,5]))
  R.init <- ceiling(exp(l2.emu.par[ii,6]))
  
  predictions <- f.emulator(x.test, l2.emu.par[ii,])
  scale.pred <- predictions$scale.pred
  pred.s[,ii] <- predictions$mean
  
  names(SEIR_par) <- c("beta", "kappa", "gamma")
  cat("stochastic: R0 =", SEIR_par[1]/SEIR_par[3], "\n")
  
  para.s[para.s$country == country,-1] <- c(I.init, E.init, R.init, SEIR_par)
  
  # compute variance of R0
  # covariance matrix of parameter estimators
  cov.mx <- var.theta(length(xp), etahat.cv, f.emulator, x.test, l2.emu.par[ii,], method="L2", emulate.fg=TRUE, side=rep(1,6))
  
  # delta method to compute standard deviation of R0
  R0.s[R0.s[,1] == country, 2] <- SEIR_par[1]/SEIR_par[3]
  R0.s[R0.s[,1] == country, 3] <- sqrt(matrix(c(1/SEIR_par[3],-SEIR_par[1]/(SEIR_par[3])^2), nrow = 1) %*% cov.mx[-c(2,4:6),-c(2,4:6)] %*% matrix(c(1/SEIR_par[3],-SEIR_par[1]/(SEIR_par[3])^2), ncol = 1))
  R0.s[R0.s[,1] == country, 3] <- R0.s[R0.s[,1] == country, 3] * sqrt(phi.hat[ii]) # times dispersion parameter
  
  # delta method to compute standard deviation of 1/kappa
  incubation[ii] <- 1/SEIR_par[2]
  incubation.sd[ii] <- sqrt(matrix(-1/(SEIR_par[2])^2, nrow = 1) %*% cov.mx[2,2] %*% matrix(-1/(SEIR_par[2])^2, ncol = 1)) * sqrt(phi.hat[ii])
  
  # delta method to compute the variance of model outputs
  logf.emulator <- function(x, theta){
    xnew <- cbind(x, matrix(rep(theta, each = length(x)), ncol = length(theta)))
    predictions <- predict(hetGP.fit, xnew)
    pred.mean <- predictions$mean
    pred.var <- predictions$sd2
    return(list(mean=pred.mean, s2=pred.var))
  }
  predictions <- logf.emulator(x.test, l2.emu.par[ii,])
  df <- jacobian(function(para) {logf.emulator(x.test, para)$mean}, l2.emu.par[ii,])
  var.s <- diag(df %*% cov.mx %*% t(df)) * phi.hat[ii]
  var.s <- var.s + predictions$s2
  # upper quantile
  uq.s[,ii] <- exp(predictions$mean-1 + sqrt(2*var.s)*erfinv(2*0.975 - 1))
  uq.s[,ii] <- uq.s[,ii] * scale.pred
  # lower quantile
  lq.s[,ii] <- exp(predictions$mean-1 + sqrt(2*var.s)*erfinv(2*0.025 - 1))
  lq.s[,ii] <- lq.s[,ii] * scale.pred
}

##### plot R0
# Figure 3: deterministic R0
pdf(file="deterministic_R0.pdf", width = 10, height = 5)
R0.d$upper <- R0.d[,2] + qnorm(0.975)*R0.d[,3]
R0.d$lower <- R0.d[,2] - qnorm(0.975)*R0.d[,3]
R0.d$rank <- rank(R0.d[,2])
p <- ggplot(R0.d, aes(x = reorder(country, R0), y = R0), show.legend = FALSE) + theme_bw(base_size=20)
p <- p + geom_point(show.legend = FALSE) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
  xlab("") + ylab(expression(R[0]))  +  
  theme(legend.position = "none",   axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off()

# Figure 5: stochastic R0
pdf(file="stochastic_R0.pdf", width = 10, height = 5)
R0.s$upper <- R0.s[,2] + qnorm(0.975)*R0.s[,3]
R0.s$lower <- R0.s[,2] - qnorm(0.975)*R0.s[,3]
R0.s$rank <- rank(R0.s[,2])
p <- ggplot(R0.s, aes(x = reorder(country, R0), y = R0), show.legend = FALSE) + theme_bw(base_size=20)
p <- p + geom_point(show.legend = FALSE) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
  xlab("") + ylab(expression(R[0]))  +  
  theme(legend.position = "none",   axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off()

# Figure 6: incubation period
tmp.df <- data.frame("country"=country.select,
                     "invkappa"=incubation,
                     "invkappa_lower"=incubation - qnorm(0.975)*incubation.sd,
                     "invkappa_upper"=incubation + qnorm(0.975)*incubation.sd)
pdf(file="stochastic_kappa.pdf", width = 10, height = 5)
p <- ggplot(tmp.df, aes(x = reorder(country, invkappa), y = invkappa), show.legend = FALSE) + theme_bw(base_size=20)
p <- p + geom_point(show.legend = FALSE) + geom_errorbar(aes(ymin = invkappa_lower, ymax = invkappa_upper), width = 0.2) + 
  xlab("") + ylab(expression(1/kappa))  +  
  theme(legend.position = "none",   axis.text.x = element_text(angle = 90, hjust = 1))
p + geom_hline(yintercept =  mean(tmp.df$invkappa), linetype="dashed", colour = "red")
dev.off()

##### plot model output results
# Figure 4: deterministic
pdf(file="deterministic_fit.pdf", width = 8, height = 6)
par(oma=c(0,0,0,0), mar=c(5.1, 2.1, 2.1, 1), mfrow=c(3,4))
for(i in 1:12){
  country <- R0.d$country[R0.d$rank == n-i+1]
  yp <- Infected.out[,country.select == country]
  xp <- 1:length(yp)
  xp <- xp[yp>0]
  yp <- yp[yp>0]
  plot(xp, yp, xaxt = "n", xlab = "", ylab = "active cases", main = country,
       col="darkgray", pch=16 , cex=1, ylim=c(0,quantile(yp,0.995)))
  axis(1, at = c(1, 32, 62, 93, 123, 154, 185, 215, 246, 276, 307, 338, 364)[seq(1,13,2)], 
       labels = c(month.name[3:12], month.name[1:3])[seq(1,13,2)], las = 2)
  lines(x.test, pred.d[,country.select == country], type = "l", lwd = 2, lty = 1, col = 2)
  lines(x.test, uq.d[,country.select == country], type = "l", lwd = 1, lty = 2, col = 2)
  lines(x.test, lq.d[,country.select == country], type = "l", lwd = 1, lty = 2, col = 2)
  polygon(c(rev(x.test), x.test), c(rev(uq.d[,country.select == country]), lq.d[,country.select == country]), col = rgb(0.9,0,0,0.1) , border = NA)
}
dev.off()

# Figure 7: stochastic
pdf(file="stochastic_fit.pdf", width = 8, height = 6)
par(oma=c(0,0,0,0), mar=c(5.1, 2.1, 2.1, 1), mfrow=c(3,4))
for(i in 1:12){
  country <- R0.s$country[R0.s$rank == n-i+1]
  yp <- Infected.out[,country.select == country]
  xp <- 1:length(yp)
  xp <- xp[yp>0]
  yp <- yp[yp>0]
  plot(xp, yp, xaxt = "n", xlab = "", ylab = "active cases", main = country,
       col="darkgray", pch=16 , cex=1, ylim=c(0,quantile(yp,0.995)))
  axis(1, at = c(1, 32, 62, 93, 123, 154, 185, 215, 246, 276, 307, 338, 364)[seq(1,13,2)], 
       labels = c(month.name[3:12], month.name[1:3])[seq(1,13,2)], las = 2)
  lines(x.test, pred.s[,country.select == country], type = "l", lwd = 2, lty = 1, col = 2)
  lines(x.test, uq.s[,country.select == country], type = "l", lwd = 1, lty = 2, col = 2)
  lines(x.test, lq.s[,country.select == country], type = "l", lwd = 1, lty = 2, col = 2)
  polygon(c(rev(x.test), x.test), c(rev(uq.s[,country.select == country]), lq.s[,country.select == country]), col = rgb(0.9,0,0,0.1) , border = NA)
}
dev.off()


