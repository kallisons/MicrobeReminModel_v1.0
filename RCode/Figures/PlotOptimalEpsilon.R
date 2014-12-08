r<-0.001
bdiam<-5E-7
blength<-1E-6
cperb<-5E-12
pom_initial<-0.6	

numB<-(4*3.14*r^2)/(bdiam*blength)
quorum<-numB*cperb

dw<-(17*(r*2*1000)^1.8)	
#pnum<-pom_initial/(dw*0.001*0.12)
pnum<-84

EOpt.ran<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Interior_LA_362to378min.txt", header=TRUE)
EOpt.int2<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Interception_LA_362to378min.txt", header=TRUE)
EOpt.ret<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Retention_LA_362to378min.txt", header=TRUE)


reduce<-as.data.frame(seq(2,110,2))
colnames(reduce)<-c("qfrac.init")
EOpt.ran<-merge(EOpt.ran, reduce, all.y=TRUE)
EOpt.int2<-merge(EOpt.int2, reduce, all.y=TRUE)

OutGraph<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Graphs/OptimalEpsilons4.ps", sep="")
postscript(OutGraph, family="Helvetica", width=3, height=3, pointsize=10)

#quartz(width=4, height=4)
par(oma=c(6.5,3,1,1))
par(mar=c(1,1,0.5,0.5)+0.1)
par(mfrow=c(1,1))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

library(RColorBrewer)
colpick<-brewer.pal(4, "Spectral")

#Epsilon Opt
#plot(EOpt.ran$qfrac.init, EOpt.ran$maxe, ylim=c(0, 0.05), xlim=c(0,110), xaxt="n", type="l", lwd=2, col="red", las=1)
plot(EOpt.ran$qfrac.init, EOpt.ran$maxe, ylim=c(0, 0.05), xlim=c(0,110), xaxt="n", type="p", col=colpick[1], pch=17, cex=1, las=1)
abline(v=25, lty=3, col="grey20")
abline(v=50, lty=3, col="grey20")
abline(v=100, lty=3, col="grey20")


axis(side=1, labels=TRUE)
per<-seq(0, 100, 20)
babund<-round(((per/100)*quorum/cperb)/1E6, digits=0)
axis(side=1, col="black", line=4, labels=babund, at=per)



##lines(EOpt.int2$qfrac.init, EOpt.int2$maxe, col="green4", lwd=2)
##lines(EOpt.ret$qfrac.init, EOpt.ret$maxe, col="blue", lwd=2)
points(EOpt.ran$qfrac.init, EOpt.ran$maxe, col=colpick[1], pch=17, cex=1)
points(EOpt.int2$qfrac.init, EOpt.int2$maxe, col=colpick[2], pch=16, cex=0.8)
lines(EOpt.ret$qfrac.init, EOpt.ret$maxe, col=colpick[4], lwd=3.5)

legend(2, 0.98*0.05, bg="white", seg.len=1, c( "Interior","Interception","Retention"), lwd=3.5, col=c("black", "black", col=colpick[4]), lty=c(2,3,1), y.intersp=0.9)

mtext(expression(paste("Optimal ",epsilon, " (", h^{-1}, ")")), side=2, line=2.75, cex=1.2)
mtext("Bacterial coverage of particle (%)", side=1, line=1.2, cex=1.2, adj=0.5, outer=TRUE)
mtext(expression(paste("Bacterial abundance  (",10^{6}, " ", particle^{-1}, ")")), side=1, line=6.4, cex=1.2, adj=0.5)



dev.off()



#abline(v=22, lwd=1)
#abline(v=32, lwd=1)
#abline(h=0.028, lwd=1)

quartz(width=4, height=4)
par(oma=c(3,3,1,1))
par(mar=c(1,1,1,1)+0.1)
par(mfrow=c(1,1))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

#plot(EOpt.ret$qfrac.init, EOpt.ret$maxe, col="green4", ylim=c(0,0.05), lwd=2, type="l")
#top<-63
#width<-75
#a<-0.074
#b<--0.017
#qfrac<-seq(1,110,1)
#ep<-a*exp(qfrac*b)*(1-((qfrac-top)/(width/2))^2)
#lines(qfrac, ep, ylim=c(0,0.05), type="l", lwd=2)
#abline(h=0.006)

#abline(v=26) #lower threshold
#abline(v=92) #upper threshold

#qfrac<-seq(1,110,1)
#emax=0.027
#et = 37
#ke = 52
#x = 1.42
#plot(EOpt.ran$qfrac.init, EOpt.ran$maxe, ylim=c(0, 0.05), xlim=c(0,110), type="l", lwd=2, col="red")
#qfrac<-seq(1,110,1)
#ep<-emax*abs(qfrac-et)^x/((ke-et)^x+abs(qfrac-et)^x)
#ep2<-c(diff(ep),1)
#ep<-ifelse(ep2<0, 0, ep)
#lines(qfrac, ep, lwd=2, col="lightgrey")


emax=0.025
et = 31
ke = 43
x = 1.44
plot(EOpt.int2$qfrac.init, EOpt.int2$maxe, ylim=c(0, 0.05), xlim=c(0,110), type="l", lwd=2, col="red")
qfrac<-seq(1,110,1)
ep<-emax*abs(qfrac-et)^x/((ke-et)^x+abs(qfrac-et)^x)
ep2<-c(diff(ep),1)
ep<-ifelse(ep2<0, 0, ep)
lines(qfrac, ep, lwd=2, col="lightgrey")




