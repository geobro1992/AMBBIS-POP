#--------------------------------------------------------------------
# Analysis of the JS model under the superpopulation parameterization
#--------------------------------------------------------------------
library(R2WinBUGS)
library(BaSTA)

#Chose Pond
pid = 4

# Import data:
dd = read.csv("raw.csv")                
dd[, 4] = as.Date(dd[, 4], format = "%m/%d/%Y")
dd = subset(dd, dd$X.10 == pid)
dd = na.omit(dd[, 3:4])
dd = dd[which(dd[, 1] != ""), ]

Y = CensusToCaptHist(ID = dd[, 1], d = dd[, 2])
CH <- as.matrix(Y[, 2:9])

# Augment capture-histories by 20, 200, & 2000 pseudo-individuals
nz <- c(20, 200, 2000)
CH.aug <- rbind(as.matrix(CH), matrix(0, ncol = dim(CH)[2], nrow = nz[1]))
CH.aug1 <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz[2]))
CH.aug2 <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz[3]))

# Bundle data
bugs.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])
bugs.data1 <- list(y = CH.aug1, n.occasions = dim(CH.aug1)[2], M = dim(CH.aug1)[1])
bugs.data2 <- list(y = CH.aug2, n.occasions = dim(CH.aug2)[2], M = dim(CH.aug2)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), beta = c(1, NA, NA, NA, 1, NA, 1, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), z = CH.aug)}  
inits1 <- function(){list(mean.phi = runif(1, 0, 1), beta = c(1, NA, NA, NA, 1, NA, 1, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), z = CH.aug1)}  
inits2 <- function(){list(mean.phi = runif(1, 0, 1), beta = c(1, NA, NA, NA, 1, NA, 1, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), z = CH.aug2)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "nu")

# MCMC settings
ni <- 500000
nt <- 100
nb <- 10000
nc <- 1

# Call WinBUGS from R (BRT 40 min)
js.super <- bugs(bugs.data, inits, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

js.super1 <- bugs(bugs.data1, inits1, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

js.super2 <- bugs(bugs.data2, inits2, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())


save(js.super, file = "POP_Output20.RData")
save(js.super1, file = "POP_Output200.RData")
save(js.super2, file = "POP_Output2000.RData")





