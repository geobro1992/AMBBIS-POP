# 10.3.3. The superpopulation parameterization
sink("js-super.bug")
cat("
model {

   for (i in 1:M) {
    for (t in 1:(n.occasions-1)) {
      phi[i, t] <- mean.phi[group[i]]
    } #t

    for (t in 1:n.occasions) {
      p[i, t] <- mean.p[group[i], t]
    } #t
  } #i

for (g in 1:G) {
  mean.phi[g] ~ dunif(0, 1)         # Prior for mean survival
  
    for (t in 1:n.occasions) {
      mean.p[g, t] ~ dunif(0, 1)           # unequal cap in 1st year
  }
}
psi ~ dunif(0, 1)              # Prior for inclusion probability
#----------------------------------------    
# Dirichlet prior for entry probabilities

beta[1] ~ dgamma(1, 1)
beta[2] ~ dgamma(1, 1)
beta[3] ~ dgamma(.01, 1)
beta[4] ~ dgamma(.01, 1)
beta[5] ~ dgamma(1, 1)
beta[6] ~ dgamma(.01, 1)
beta[7] ~ dgamma(.01, 1)
beta[8] ~ dgamma(1, 1)

  for (t in 1:n.occasions) {
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
  group[i] ~ dcat(pp[]) 
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

pp[1] <- .6
pp[2] <- .4
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

for (i in 1:M) {
  for (t in 1:n.occasions) {

    N1[i, t] <- u[i, t] * (1-(group[i]-1))        # Actual population size
    B1[i, t] <- recruit[i, t] * (1-(group[i]-1))  # Number of entries

    N2[i, t] <- u[i, t] * (group[i]-1)        # Actual population size
    B2[i, t] <- recruit[i, t] * (group[i]-1)  # Number of entries

  } #t
} #i

for (t in 1:n.occasions) {
 
 N4[t] <- sum(N1[1:M, t]) 
 N5[t] <- sum(N2[1:M, t]) 

 B4[t] <- sum(B1[1:M, t]) 
 B5[t] <- sum(B2[1:M, t]) 

}

  for (i in 1:M) {
    Nind[i] <- sum(u[i, 1:n.occasions])
    Nalive[i] <- 1 - equals(Nind[i], 0)
  } #i
  
  Nsuper <- sum(Nalive[])         # Superpopulation size

}
  ", fill = TRUE)
sink()
