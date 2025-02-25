load(file="res.rho.Diet.RData")
load(file="res.rho.ORIGIN.RData")
load(file="res.rho.sim001.n152.RData")
res.Com2seq.CS.sim001 <- res.Com2seq.CS.sim01 # correct typo
load(file="res.rho.sim01.n76.RData")


pdf("res_rho.pdf", height=5, width=6)
par(mfrow=c(1,1), mar=c(4.5,5,4,0.5))

plot(density(res.Com2seq.CS.ORIGIN$rho_estimated_all), ylim=c(0, 3), col="blue", main="16S-SMS data correlations", xlab=expression(hat(rho)), ylab="Density", lwd=2)
lines(density(res.Com2seq.CS.sim001$rho_estimated_all), lty=2, col="blue", lwd=2)
lines(density(res.Com2seq.CS.Diet$rho_estimated_all), col="red", lwd=2)
lines(density(res.Com2seq.CS.sim01$rho_estimated_all), lty=2, col="red", lwd=2)
legend(x="topright", legend=c("ORIGIN data", expression(paste("Sim, ", tau, " = 0.001")), "Diet data", expression(paste("Sim, ", tau, " = 0.01"))), 
       col=c("blue", "blue", "red", "red"), lty=c(1,2,1,2), bty="n", lwd=c(2,2,2,2))

dev.off()










