SEIR.sim <- function(n.rep, n.train, N, xp, lower, upper){
  # n.rep: number of replicates for training the stochastic SEIR model
  # n.train: sample size of training data for the stochastic SEIR model
  
  xs.train <- unique(c(seq(1,length(xp),20), length(xp)))
  n.timestep <- length(xs.train)
  
  set.seed(1)
  d <- 6
  theta.train <- maximinLHS(n = n.train, k = d)
  theta.train <- t(t(theta.train) * (upper - lower) + lower)
  
  sim.out <- matrix(0, nrow = n.train * n.rep, ncol = n.timestep)
  start.time <- Sys.time()
  for(jj in 1:n.train){
    print(jj)
    beta <- theta.train[jj,1]
    epsilon <- theta.train[jj,2]
    gamma <- theta.train[jj,3]
    I.init <- ceiling(exp(theta.train[jj,4]))
    E.init <- ceiling(exp(theta.train[jj,5]))
    R.init <- ceiling(exp(theta.train[jj,6]))
    
    model <- SEIR(u0 = data.frame(S = rep(N - I.init - E.init - R.init, n.rep), 
                                  E = rep(E.init, n.rep), I = rep(I.init, n.rep), R = rep(R.init, n.rep)),
                  tspan = xs.train, beta = beta, epsilon = epsilon, gamma = gamma)
    result <- run(model)
    result.out <- trajectory(result)[,c(4,2)]
    result.out[,1] <- result.out[,1] * epsilon
    out <- unstack(result.out)
    sim.out[((jj-1)*n.rep+1) :(jj*n.rep),] <- as.matrix(out)
  }
  end.time <- Sys.time()
  time.proc <- difftime(end.time, start.time)
  
  Xs <- matrix(rep(theta.train, each = n.rep * n.timestep), ncol = ncol(theta.train))
  Xs <- cbind(rep(xs.train, n.rep * n.train), Xs)
  Ys <- c(t(sim.out))
  
  return(list(Xs=Xs, Ys=Ys, time.proc=time.proc))
}


