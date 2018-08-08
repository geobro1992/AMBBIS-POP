#--------------------------------------------------------------------
# Analysis of the JS model under the superpopulation parameterization
#--------------------------------------------------------------------
library(R2WinBUGS)
library(BaSTA)

# DATA PREP.:
# Import data:
r.data = read.csv("raw.csv")                
r.data[, 4] = as.Date(r.data[, 4], format = "%m/%d/%Y")

# Assign codes to indv without Master IDs
r.data = na.omit(r.data[, c(3:4, 12, 19)])
r.data[, 1] = as.character(r.data[, 1])

count = 2000

for (i in 1:length(r.data[,1])) {
  if (r.data[i, 1] == "") {
    r.data[i, 1] = count
    count = count + 1
  }
}

colnames(r.data) = c("ID", "Date", "Pond","Age")
r.data[,1] = as.factor(r.data[,1])

r.data = r.data[which(r.data$Age != "M"),]
r.data = r.data[which(r.data$Pond != "53"),]

Y = CensusToCaptHist(ID = r.data[, 1], d = r.data[, 2])

# get pond ids for group effect
pid = unique(r.data[,c(1,3)])
n_occur <- data.frame(table(pid[,1]))
n_occur[n_occur$Freq > 1,] # find duplicates

pid = pid[-c(16, 291, 434, 539, 632), ]

pid[, 1] = as.numeric(pid[, 1])
pid[, 2] = as.factor(as.numeric(pid[, 2]))
pid = data.frame(pid[,2])
colnames(pid) = "Pond"

# remove pond 53 individuals
pid = append(pid[,1], rep(NA, 2000))
# Augment capture-histories by 20, 200, & 2000 pseudo-individuals
nz <- c(20, 200, 2000)

CH <- as.matrix(Y[, 2:9])

CH.aug <- rbind(as.matrix(CH), matrix(0, ncol = dim(CH)[2], nrow = nz[1]))
CH.aug1 <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz[2]))
CH.aug2 <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz[3]))

# Bundle data
bugs.data <- list(y = CH, n.occasions = dim(CH)[2], M = dim(CH)[1], G = 2, group = pid)
bugs.data1 <- list(y = CH.aug1, n.occasions = dim(CH.aug1)[2], M = dim(CH.aug1)[1], G = 2, group = pid)
bugs.data2 <- list(y = CH.aug2, n.occasions = dim(CH.aug2)[2], M = dim(CH.aug2)[1], G = 2, group = pid)

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), beta = c(1, NA, NA, NA, 1, NA, 1, 1), mean.p = runif(2, 0.5, 1), psi = runif(1, 0, 1), z = CH)}  
inits1 <- function(){list(mean.phi = runif(2, 0, 1), beta = c(1, NA, NA, NA, 1, NA, 1, 1), mean.p = runif(2, 0.5, 1), psi = runif(1, 0, 1), z = CH.aug1)}  
inits2 <- function(){list(mean.phi = runif(2, 0, 1), beta = c(1, 1, 1, 1, 1, 1, 1, 1), group = c(rep(NA, 716), rep(1, 2000)),
                          mean.p = rbind(runif(8, 0.5, 1), runif(8, 0.5, 1)), psi = runif(1, 0, 1), z = CH.aug2)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "Nsuper", "B4", "B5", "N4", "N5")

# MCMC settings
ni <- 100000
nt <- 90
nb <- 10000
nc <- 3

# Call WinBUGS from R
js.super <- bugs(bugs.data, inits, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

js.super1 <- bugs(bugs.data1, inits1, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

js.super2 <- bugs(bugs.data2, inits2, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                 debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd(), DIC = F)


save(js.super, file = "POP_Output20.RData")
save(js.super1, file = "POP_Output200.RData")
save(js.super2, file = "POP_Output_Jul_20_18.RData")





