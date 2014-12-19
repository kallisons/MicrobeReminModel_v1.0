#Change working directory to MicrobeReminModel_v1.0 using the setwd() command or using the Misc menu.  
getwd() #shows the current working directory for the R gui.

#----------------------
# Read in Model Output
#----------------------

bfolder<-paste("../Output/MRM1D/BacteriaRates/")
bfiles<-list.files(bfolder)

for(i in 1:length(bfiles)){

brates<-read.table(paste(bfolder, bfiles[i], sep=""))
colnames(brates)<-c("depth", "attach", "detach", "ba_grow", "bf_grow", "b_ratio", "bp_ratio", "up_ba", "up_bf")

brates$depth<-abs(brates$depth)*-1

#----------------
#  Create Plot
#----------------

#OutGraph<-paste("/Users/kasmith/Dropbox/OCB2012/Poster/Rates.ps", sep="")
#postscript(OutGraph, family="Helvetica", width=7.5, height=6)

#quartz(width=7.5, height=6)
par(oma=c(2,2,2,2))
par(mar=c(1,4,4,1))
par(mfrow=c(2,3))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

pdepth<--4000

#--------------------
# Plot Bacteria Rates
#--------------------

#The rates are multiplied by 24 to transform from hourly to daily rates

plot(brates$attach*24, brates$depth, type = "l", ylim=c(pdepth, 0), xlim=c(0,(max(brates$attach*24)+max(brates$attach)*24*0.10)), lwd=2, xaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("Attachment (", d^-1, ")", sep="")), line=2.7)
axis(side=3)

plot(brates$up_bf*24, brates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(brates$up_bf)*24+max(brates$up_bf)*24*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Uptake DB (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(brates$bf_grow*24, brates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(brates$bf_grow)*24+max(brates$bf_grow)*24*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Growth DB (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(brates$detach*24, brates$depth, type = "l", ylim=c(pdepth, 0), xlim=c(0,(max(brates$detach*24)+max(brates$detach)*24*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Detachment (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(brates$up_ba*24, brates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(brates$up_ba)*24+max(brates$up_ba)*24*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Uptake PB (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(brates$ba_grow*24, brates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(brates$ba_grow)*24+max(brates$ba_grow)*24*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Growth PB (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

dev.off()

}