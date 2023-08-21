### SEIR model 
SEIR.ode <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dE <- beta/N * I * S - kappa * E
    dI <- kappa * E - gamma * I
    dR <- gamma * I
    list(c(dS, dE, dI, dR))
  })
}