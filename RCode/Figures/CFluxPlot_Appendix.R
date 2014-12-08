state3<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/CFlux1D/Retention/StateVariables/StateVariables_Sink2.00_Size0.0010_Qfrac1.00.out", sep="")

#---------------------------------------
#	RETENTION: PLOT FLUXES SEPARATELY
#----------------------------------------

colnames(state3)<-c("depth", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom", "pomflux", "baflux", "enzflux", "hflux", "denzflux")
state3$flux<-state3$pomflux+state3$baflux
state3$flux<-state3$flux*24
state3$flux2<-state3$pomflux+state3$baflux+state3$enzflux+state3$hflux+state3$denzflux
state3$flux2<-state3$flux2*24
state3$depth<-abs(state3$depth)*-1
state3$pomflux<-state3$pomflux*24
state3$baflux<-state3$baflux*24
state3$enzflux<-state3$enzflux*24
state3$hflux<-state3$hflux*24
state3$denzflux<-state3$denzflux*24



#print(sum((state3$flux-min(state3$flux))*state3$depth)/sum(state3$flux-min(state3$flux)))
#print(sum((state3$flux2-min(state3$flux2))*state3$depth)/sum(state3$flux2-min(state3$flux2)))


#------------------------------
#   MAKE PLOT
#------------------------------
OutGraph<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Graphs/Figure_SeparateFluxes2.ps", sep="")
postscript(OutGraph, family="Helvetica", width=6.5, height=3, pointsize=12)

#quartz(width=6,height=3)
par(mfrow=c(1,5))
par(oma=c(3,3.5,3,1))
par(mar=c(0.5,1.5,1.5,0)+0.1)
par(xaxs="i")
par(yaxs="i")

#library(RColorBrewer)
#colpick<-brewer.pal(4, "Spectral")
#colpick<-brewer.pal(4, "RdYlBu")
#colpick<-brewer.pal(4, "Dark2")

plot(state3$pomflux, abs(state3$depth)*-1, ylim=c(-1000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col="red", type="l", lwd=2, cex=1)
lines(state3$flux, abs(state3$depth)*-1, lty=1, lwd=2, col="grey60")
lines(state3$flux2, abs(state3$depth)*-1, lty=3, lwd=2, col="grey60")
lines(state3$pomflux, abs(state3$depth)*-1, lty=1, lwd=2, col="red")
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0, 10, 20, 30), at=c(0, 10, 20, 30))
axis(side=2, las=2, tick=TRUE, labels=c(0,250,500,750,1000), at=c(0,-250,-500,-750,-1000), cex.axis=1.2)
mtext("Depth (m)", side=2, line=3.4, cex=1)
mtext("Particle", side=1, line=1, cex=1)

plot(state3$baflux, abs(state3$depth)*-1, ylim=c(-1000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col="red", type="l", lwd=2, cex=1)
lines(state3$flux, abs(state3$depth)*-1, lty=1, lwd=2, col="grey60")
lines(state3$flux2, abs(state3$depth)*-1, lty=3, lwd=2, col="grey60")
lines(state3$baflux, abs(state3$depth)*-1, lty=1, lwd=2, col="red")
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0, 10, 20, 30), at=c(0, 10, 20, 30))
axis(side=2, las=2, tick=TRUE, labels=FALSE, at=c(0,-250,-500,-750,-1000))
mtext("Attached", side=1, line=1, cex=1)
mtext("Bacteria", side=1, line=2.4, cex=1)

plot(state3$hflux, abs(state3$depth)*-1, ylim=c(-1000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col="red", type="l", lwd=2, cex=1)
lines(state3$flux, abs(state3$depth)*-1, lty=1, lwd=2, col="grey60")
lines(state3$flux2, abs(state3$depth)*-1, lty=3, lwd=2, col="grey60")
lines(state3$hflux, abs(state3$depth)*-1, lty=1, lwd=2, col="red")
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0, 10, 20, 30), at=c(0, 10, 20, 30))
axis(side=2, las=2, tick=TRUE, labels=FALSE, at=c(0,-250,-500,-750,-1000))
mtext("Hydrolysate", side=1, line=1, cex=1)

plot(state3$enzflux, abs(state3$depth)*-1, ylim=c(-1000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col="red", type="l", lwd=2, cex=1)
lines(state3$flux, abs(state3$depth)*-1, lty=1, lwd=2, col="grey60")
lines(state3$flux2, abs(state3$depth)*-1, lty=3, lwd=2, col="grey60")
lines(state3$enzflux, abs(state3$depth)*-1, lty=1, lwd=2, col="red")
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0, 10, 20, 30), at=c(0, 10, 20, 30))
axis(side=2, las=2, tick=TRUE, labels=FALSE, at=c(0,-250,-500,-750,-1000))
mtext("Active", side=1, line=1, cex=1)
mtext("Exoenzyme", side=1, line=2.4, cex=1)

plot(state3$denzflux, abs(state3$depth)*-1, ylim=c(-1000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col="red", type="l", lwd=2, cex=1)
lines(state3$flux, abs(state3$depth)*-1, lty=1, lwd=2, col="grey60")
lines(state3$flux2, abs(state3$depth)*-1, lty=3, lwd=2, col="grey60")
lines(state3$denzflux, abs(state3$depth)*-1, lty=1, lwd=2, col="red")
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0, 10, 20, 30), at=c(0, 10, 20, 30))
axis(side=2, las=2, tick=TRUE, labels=FALSE, at=c(0,-250,-500,-750,-1000))
mtext("Inactive", side=1, line=1, cex=1)
mtext("Exoenzyme", side=1, line=2.4, cex=1)

xaxis<-expression(paste("Carbon flux (mg ", m^-2," ", d^-1, ")", sep=""))
mtext(xaxis, side=3, line=1, cex=1, outer=TRUE)

dev.off()


