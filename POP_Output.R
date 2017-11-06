
load(file = "POP_Output20.RData")
POP_Output20 = js.super

load(file = "POP_Output200.RData")
POP_Output200 = js.super1

load(file = "POP4_Output2000.RData")
POP_Output2000 = js.super2 

print(POP_Output20, digits = 2)
print(POP_Output200, digits = 2)
print(POP_Output2000, digits = 2)

# Code to produce Fig
par(mfrow = c(1,1), mar = c(5, 6, 2, 1), mgp = c(3.4, 1, 0), las = 1)
plot(density(POP_Output20$sims.list$Nsuper), main = "", xlab = "", ylab = "Density", frame = FALSE, lwd = 2, ylim = c(0, 0.023), col = "blue")
abline(v = POP_Output20$mean$Nsuper, col = "red", lwd = 2)
mtext("Size of Superpopulation", 1, line = 3)

time <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017)

b3.lower <- b3.upper <- numeric()
for (t in 1:4) {
  b3.lower[t] <- quantile(POP_Output20$sims.list$b[,t], 0.025)
  b3.upper[t] <- quantile(POP_Output20$sims.list$b[,t], 0.975)
}

plot(x = time + .5, y = POP_Output20$mean$b,
     main = "Pond 4", xlab = "", ylab = "Entry probability",
     frame = FALSE, las = 1, pch = 16,
     xlim = c(2010, 2018), ylim = c(0, max(b3.upper) + 0.1))
segments(time + .5, b3.lower, time + .5, b3.upper)
mtext("Year", 1, line = 3)


N.lower <- N.upper <- numeric()
for (t in 1:length(time)) {
  N.lower[t] <- quantile(POP_Output200$sims.list$N[,t], 0.025)
  N.upper[t] <- quantile(POP_Output200$sims.list$N[,t], 0.975)
}

plot(x = time + .5, y = POP_Output200$mean$N, xlab = "", ylab = "Population Size", main = "Pond 5",
     frame = FALSE, las = 1, pch = 16,
     xlim = c(2010, 2018), ylim = c(0, max(N.upper) + 100))
segments(time + .5, N.lower, time + .5, N.upper)
mtext("Year", 1, line = 3)

POP_Output20$mean$mean.p
POP_Output20$mean$mean.phi

# model validation - aug sensitivity 
N1.lower <- N2.lower <- N3.lower <- N1.upper <- N2.upper <- N3.upper <- numeric()
for (t in 1:length(time)) {
  N1.lower[t] <- quantile(POP_Output20$sims.list$N[,t], 0.025)
  N2.lower[t] <- quantile(POP_Output200$sims.list$N[,t], 0.025)
  N3.lower[t] <- quantile(POP_Output2000$sims.list$N[,t], 0.025)
  N1.upper[t] <- quantile(POP_Output20$sims.list$N[,t], 0.975)
  N2.upper[t] <- quantile(POP_Output200$sims.list$N[,t], 0.975)
  N3.upper[t] <- quantile(POP_Output2000$sims.list$N[,t], 0.975)
}

time <- time + 0.5

pdf("P4_Pop_Plot17.pdf", 12, 8)
plot(x = time - 0.15, y = POP_Output20$mean$N, xlab = "Year", ylab = "Population Size", 
     frame = FALSE, las = 1, pch = 16, xlim = c(2010, 2018),  ylim = c(0, 1000), cex.axis = 1.3, cex.lab = 1.5)
segments(time - 0.15, N1.lower, time - 0.15, N1.upper)
points(x = time, y = POP_Output200$mean$N, pch = 1)
segments(time, N2.lower, time, N2.upper)
points(x = time + 0.15, y = POP_Output2000$mean$N, pch = 17)
segments(time + 0.15, N3.lower, time + 0.15, N3.upper)
legend("topright", legend = c("+20", "+200", "+2000"), pch = c(16, 1, 17))
dev.off()


pdf("P4_Bugs_Plot+20.pdf", 12, 15)
plot(POP_Output20)
dev.off()

pdf("P4_Bugs_Plot+200.pdf", 12, 15)
plot(POP_Output200)
dev.off()

pdf("P4_Bugs_Plot+2000.pdf", 12, 15)
plot(POP_Output2000)
dev.off()



# Diagnostics to try
library(coda) 

outmc <- as.mcmc.list(POP_Output20)
summary(outmc)
plot(outmc, ask=TRUE)
outmc[,"alpha"]
autocorr.plot(outmc, ask=TRUE)
gelman.plot(outmc, ask=TRUE)
gelman.diag(outmc)
window(outmc, thin=5) 

