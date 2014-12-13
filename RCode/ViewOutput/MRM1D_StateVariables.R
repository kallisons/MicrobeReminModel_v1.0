#Change working directory to MicrobeReminModel_v1.0 using the setwd() command or using the Misc menu.  
getwd() #shows the current working directory for the R gui.  

#----------------------
# Read in Model Output
#----------------------

state<-read.table("Output/MRM1D/StateVariables/StateVariables_Sink2.00_Size0.0010_Qfrac0.25.out")

colnames(state)<-c("depth", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom", "pomflux", "baflux")

state$depth<-state$depth*-1

OutGraph<-paste("Graphs/MRM1D1_StateVariables.ps", sep="")
#Note preserve width/height aspect ratio at 0.8
#postscript(OutGraph, family="Helvetica", width=7.2, height=6)
#postscript(OutGraph, family="Helvetica", width=6, height=5)

#----------------
# Create Plot
#----------------

quartz(width=6, height=6.5)
par(oma=c(2,2,2,2))
par(mar=c(1,4,4,1)+0.1)
par(mfrow=c(3,3))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

pdepth<--4000

#----------------------
# Plot State Variables
#----------------------

plot(state$pom, state$depth, type = "l", ylim=c(pdepth, 0), xlim=c(0,(max(state$pom)+max(state$pom)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("POM (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, las=2, tick=TRUE, labels=seq(0, -1*pdepth, 500), at=seq(0,pdepth,-500))

plot(state$enz, state$depth, type = "l", ylim=c(pdepth, 0), xlim=c(0,(max(state$enz)+max(state$enz)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Enzyme Part. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$ba, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$ba)+max(state$ba)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Attached Bact. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$hdom, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$hdom)+max(state$hdom)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("Hydrolysate Diss. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, las=2, tick=TRUE, labels=seq(0, -1*pdepth, 500), at=seq(0,pdepth,-500))


plot(state$h, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$h)+max(state$h)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Hydrolysate Part. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$bf, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$bf)+max(state$bf)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Free-living Bact. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$denz, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$denz)+max(state$denz)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Dead Enz. Part. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$edom, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$edom)+max(state$edom)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Active Enz. Diss. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(state$dedom, state$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(state$dedom)+max(state$dedom)*0.10)), lwd=2, col="black", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Dead Enz. Diss. (mg C ", m^-3, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

#dev.off()

