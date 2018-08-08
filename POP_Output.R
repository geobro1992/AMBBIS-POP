library(boot)

load(file = "POP_Output20.RData")
POP_Output20 = js.super

load(file = "POP_Output200.RData")
POP_Output200 = js.super1

load(file = "POP_Output_Apr_25_18.RData")
POP_Output2000 = js.super2 

print(POP_Output20, digits = 2)
print(POP_Output200, digits = 2)
print(POP_Output2000, digits = 2)


# check convergence
outmc = as.mcmc.list(POP_Output2000)
summary(outmc)
plot(outmc, ask = T)
outmc[,"alpha"]
autocorr.plot(outmc, ask = T)
gelman.plot(outmc, ask = T)
gelman.diag(outmc)
window(outmc, thin = 5)


# Code to produce Fig
par(mfrow = c(1,1), mar = c(5, 6, 2, 1), mgp = c(3.4, 1, 0), las = 1)
plot(density(POP_Output2000$sims.list$Nsuper), main = "", xlab = "", ylab = "Density", frame = FALSE, lwd = 2, ylim = c(0, 0.023), col = "blue")
abline(v = POP_Output2000$mean$Nsuper, col = "red", lwd = 2)
mtext("Size of Superpopulation", 1, line = 3)

time <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017)

b4.lower <- b4.upper <- b5.lower <- b5.upper <- numeric()
for (t in 1:8) {
  b4.lower[t] <- quantile(POP_Output2000$sims.list$B4[,t], 0.025)
  b4.upper[t] <- quantile(POP_Output2000$sims.list$B4[,t], 0.975)
  b5.lower[t] <- quantile(POP_Output2000$sims.list$B5[,t], 0.025)
  b5.upper[t] <- quantile(POP_Output2000$sims.list$B5[,t], 0.975)
}

par(mfrow = c(2,2))
plot(x = time + .5, y = POP_Output2000$mean$B4,
     main = "Pond 4", xlab = "", ylab = "Entries",
     frame = FALSE, las = 1, pch = 16,
     xlim = c(2010, 2018), ylim = c(0, max(b4.upper) + 0.1))
segments(time + .5, b4.lower, time + .5, b4.upper)
mtext("Year", 1, line = 3)

plot(x = time + .5, y = POP_Output2000$mean$B5,
     main = "Pond 5", xlab = "", ylab = "Entries",
     frame = FALSE, las = 1, pch = 16,
     xlim = c(2010, 2018), ylim = c(0, max(b5.upper) + 0.1))
segments(time + .5, b5.lower, time + .5, b5.upper)
mtext("Year", 1, line = 3)


N4.lower <- N4.upper <- N5.lower <- N5.upper <- numeric()
for (t in 1:length(time)) {
  N4.lower[t] <- quantile(POP_Output2000$sims.list$N4[,t], 0.025)
  N4.upper[t] <- quantile(POP_Output2000$sims.list$N4[,t], 0.975)
  N5.lower[t] <- quantile(POP_Output2000$sims.list$N5[,t], 0.025)
  N5.upper[t] <- quantile(POP_Output2000$sims.list$N5[,t], 0.975)
  
  }

par(mfrow = c(1,1), mar = c(6, 6, 3, 3))
plot(x = time + .45, y = POP_Output2000$mean$N4, xlab = "", ylab = "", main = "",
     frame = FALSE, las = 1, pch = 16, axes = F,
     xlim = c(2010, 2018), ylim = c(0, max(N4.upper) + 100))
segments(time + .45, N4.lower, time + .45, N4.upper)

points(x = time + .55, y = POP_Output2000$mean$N5, las = 1, pch = 17, col = "dark grey")
segments(time + .55, N5.lower, time + .55, N5.upper, col = "dark grey")

axis(1, at = time + .5, labels = c("10-11", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-18"))
axis(2, at = c(0, 100, 200, 300, 400, 500, 600, 700))

mtext("Breeding Season", 1, line = 3, cex = 2)
mtext("Population Size", side = 2, las = 3, line = 3, cex = 2)

legend(2010.5, 600, 
       legend = c(expression("1; N"[e] == 204 (176-239)), 
                  expression("2; N"[e] == 130 (100-165))), 
                  pch = c(16, 17), 
       col = c("black", "dark grey"), box.lty = 0, cex = 1.5)

#------------------------------
# Effective population sizes
#------------------------------
ne4 = 1 / (sum(1/POP_Output2000$mean$N4) / length(time))
ne4u = 1 / (sum(1/N4.upper) / length(time))
ne4l = 1 / (sum(1/N4.lower) / length(time))

ne5 = 1 / (sum(1/POP_Output2000$mean$N5) / length(time))
ne5u = 1 / (sum(1/N5.upper) / length(time))
ne5l = 1 / (sum(1/N5.lower) / length(time))
#------------------------------
# Synchrony tests

cor(POP_Output2000$mean$N4, POP_Output2000$mean$N5, method = "pearson")
cor.test(POP_Output2000$mean$N4, POP_Output2000$mean$N5, method = "pearson")

cor(POP_Output2000$mean$N4, POP_Output2000$mean$N5, method = "kendall")
cor.test(POP_Output2000$mean$N4, POP_Output2000$mean$N5, method = "kendall")

df = cbind(mu.x = POP_Output2000$mean$N4, mu.y = POP_Output2000$mean$N5)

x=vector()
y=vector()

cor.mu.sigma <- function(df, n) {
     df = df[n,]
      for(i in 1:8){
          x[i] <- sample(POP_Output2000$sims.list$N4[,i], size = 1)
          y[i] <- sample(POP_Output2000$sims.list$N5[,i], size = 1)
      }
     return(cor(x, y, method = "kendall"))
}

boot.cor.mu.sigma = boot(df, cor.mu.sigma, R = 100000)
boot.ci(boot.cor.mu.sigma, type = c("norm", "basic", "perc"))

#----------------------------------
POP_Output2000$mean$mean.p
POP_Output2000$mean$mean.phi


P4 = which(POP_Output2000$mean$p[,1] > 0.5)
P5 = which(POP_Output2000$mean$p[,1] < 0.5)


x = POP_Output2000$sims.list$p[,P4,]
dim(x) <- c(1000 * 1062 , 8)

boxplot(x)

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

outmc <- as.mcmc.list(POP_Output2000)
summary(outmc)
plot(outmc, ask=TRUE)
outmc[,"alpha"]
autocorr.plot(outmc, ask=TRUE)
gelman.plot(outmc, ask=TRUE)
gelman.diag(outmc)
window(outmc, thin=5) 

# Goodness-of-Fit assessment in Bayesian analyses

par(mfrow=c(1, 2))
plot(outmc$mean$predicted, outmc$mean$residual, main="Residuals vs. predicted values", las=1, xlab="Predicted values", ylab="Residuals", pch=20, col="blue")
abline(h=0)

plot(x, out$mean$residual, main="Residuals vs. predictor values", las=1, xlab="Predictor values", ylab="Residuals", pch=20, col="blue")
abline(h=0)
par(mfrow=c(1, 1))


# Assessment of model adequacy based on posterior predictive distributions:
# 1. graphically, in a plot of the lack of fit for the ideal data vs. the lack of fit for the actual data 
## -- if the model fits the data, then about half of the points should lie above and half of them below a 1:1 line

lim <- c(0, max(c(out$sims.list$fit, out$sims.list$fit.new)))
plot(out$sims.list$fit, out$sims.list$fit.new, main="Graphical posterior predictive check", las=1, xlab="SSQ for actual data set", ylab="SSQ for ideal (new) data sets", xlim=lim, ylim=lim, col="blue")
abline(0, 1)

# 2. by means of a numerical summary called a Bayesian p-value
## -- this quantifies the proportion of times when the discrepancy measure for the perfect data sets is greater than the discrepancy measure computed for the actual data set
## -- a fitting model has a Bayesian p-value near 0.5 and values close to 0 or close to 1 suggest doubtful fit of the model

mean(out$sims.list$fit.new>out$sims.list$fit)


plot(1:10, (1:10)^2, xaxt = 'n', yaxt = 'n', 
     xlab = "Isolation", ylab = "Extinction Rate",
     type = "l")
lines(1:10, (10:1)^2, lty = "dashed")
text(2.5, 90, "Global")
text(8.5, 90, "Local")
