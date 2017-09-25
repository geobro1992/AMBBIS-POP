# 10.3.3. The superpopulation parameterization
sink("js-super.bug")
cat("
model {
    
  for (i in 1:M) {
    for (t in 1:(n.occasions-1)) {
      phi[i, t] <- mean.phi
    } #t
    
    for (t in 1:n.occasions) {
      p[i, t] <- mean.p
    } #t
  } #i
    
  mean.phi ~ dunif(0, 1)         # Prior for mean survival
  mean.p ~ dunif(0, 1)           # Prior for mean capture
  psi ~ dunif(0, 1)              # Prior for inclusion probability

#----------------------------------------    
# Dirichlet prior for entry probabilities
  
  for (t in 1:n.occasions) {
    beta[t] ~ dgamma(1, 1)
    b[t] <- beta[t] / sum(beta[1:n.occasions])
  }
    
#-----------------------------------------------
# Convert entry probs to conditional entry probs
  nu[1] <- b[1]
  
  for (t in 2:n.occasions){
    nu[t] <- b[t] / (1 - sum(b[1:(t-1)]))
  } #t
    
#-----------------------------------------------
# Likelihood function
  
  for (i in 1:M) { # First occasion
# State process
    w[i] ~ dbern(psi)
    z[i, 1] ~ dbern(nu[1])
# Observation process
    mu1[i] <- z[i, 1] * p[i, 1] * w[i]
    y[i, 1] ~ dbern(mu1[i])
    
    for (t in 2:n.occasions) { # Subsequent occasions
# State process
      q[i, t-1] <- 1 - z[i, t-1]
      mu2[i, t] <- phi[i, t-1] * z[i, t-1] + nu[t] * prod(q[i, 1:(t-1)]) 
      z[i, t] ~ dbern(mu2[i, t])
# Observation process
      mu3[i, t] <- z[i, t] * p[i, t] * w[i]
      y[i, t] ~ dbern(mu3[i, t])
    } #t
  } #i 
    
#----------------------------------------
# Calculate derived population parameters

  for (i in 1:M) {
    for (t in 1:n.occasions) {
      u[i, t] <- z[i, t] * w[i]     # Deflated latent state (u)
    }
  }
    
  for (i in 1:M) {                                 # First occasion
    recruit[i, 1] <- u[i, 1]
    for (t in 2:n.occasions) {                     # Subsequent Occasions
    recruit[i, t] <- (1 - u[i, t-1]) * u[i, t]
    } #t
  } #i

  for (t in 1:n.occasions) {
    N[t] <- sum(u[1:M, t])        # Actual population size
    B[t] <- sum(recruit[1:M, t])  # Number of entries
  } #t
    
  for (i in 1:M) {
    Nind[i] <- sum(u[i, 1:n.occasions])
    Nalive[i] <- 1 - equals(Nind[i], 0)
  } #i
  
  Nsuper <- sum(Nalive[])         # Superpopulation size
  
}
  ", fill = TRUE)
sink()