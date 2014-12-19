#Change working directory to MicrobeReminModel_v1.0 using the setwd() command or using the Misc menu.  
getwd() #shows the current working directory for the R gui.

#----------------------
# Read in Model Output
#----------------------

bfolder<-paste("../Output/MRM1D/BacteriaRates/")
bfiles<-list.files(bfolder)

mtfolder<-paste("../Output/MRM1D/MassTransferRates/")
mtfiles<-list.files(mtfolder)

for(i in 1:length(bfiles)){

brates<-read.table(paste(bfolder, bfiles[i]), sep="")


colnames(brates)<-c("depth", "attach", "detach", "ba_grow", "bf_grow", "b_ratio", "bp_ratio", "up_ba", "up_bf")
colnames(mtrates)<-c("depth", "pnum", "pom_dgd", "h_dom", "ehl_k", "enz_dom", "denz_dom", "bacover") 
brates$depth<-brates$depth*-1
mtrates$depth<-mtrates$depth*-1

#----------------
#  Create Plot
#----------------

#OutGraph<-paste("/Users/kasmith/Dropbox/OCB2012/Poster/Rates.ps", sep="")
#postscript(OutGraph, family="Helvetica", width=14, height=10)

quartz(width=7.5, height=6)
par(oma=c(2,2,2,2))
par(mar=c(1,4,4,1))
par(mfrow=c(2,4))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

#--------------------
# Plot Bacteria Rates
#--------------------

plot(brates$attach*24, brates$depth, type = "l", ylim=c(-2000, 0), xlim=c(0,(max(brates$attach*24)+max(brates$attach)*24*0.10)), lwd=2, col="red", xaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("Attachment (", d^-1, ")", sep="")), line=2.7)
axis(side=3)

plot(brates$detach*24, brates$depth, type = "l", ylim=c(-2000, 0), xlim=c(0,(max(brates$detach*24)+max(brates$detach)*24*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Detachment (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(brates$ba_grow*24, brates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(brates$ba_grow)*24+max(brates$ba_grow)*24*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Growth BA (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)


plot(brates$bf_grow*24, brates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(brates$bf_grow)*24+max(brates$bf_grow)*24*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Growth BF (", d^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

#plot(brates$b_ratio, brates$depth, type = "l", ylim=c(-3000, 0),  xlim=c(0,(max(brates$b_ratio)+max(brates$b_ratio)*0.10)), lwd=2, col="red")
#
#plot(brates$bp_ratio, brates$depth, type = "l", ylim=c(-3000, 0),  xlim=c(0,(max(brates$bp_ratio)+max(brates$bp_ratio)*0.10)), lwd=2, col="red")

#plot(mtrates$pnum, mtrates$depth, type = "l", ylim=c(-3000, 0),  xlim=c(0,(max(mtrates$pnum)+max(mtrates$pnum)*0.10)), lwd=2, col="red")
#
#
#plottitle<-expression(paste("Bacterial Rates (mg ", m^-3," ", h^-1, ")", sep=""))
#title(plottitle, cex.main=1.5, line=0.2, outer=TRUE)

#-------------------------
# Plot Mass Transfer Rates
#--------------------------

plot(mtrates$pom_dgd, mtrates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(mtrates$pom_dgd)+max(mtrates$pom_dgd)*0.10)), lwd=2, col="red", xaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("Hydrolysis (", h^-1, ")", sep="")), line=2.7)
axis(side=3)

plot(mtrates$h_dom, mtrates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(mtrates$h_dom)+max(mtrates$h_dom)*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Mass Transfer H (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(mtrates$enz_dom, mtrates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(mtrates$enz_dom)+max(mtrates$enz_dom)*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Mass Transfer EZ (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(mtrates$ehl_k, mtrates$depth, type = "l", ylim=c(-2000, 0),  xlim=c(0,(max(mtrates$ehl_k)+max(mtrates$ehl_k)*0.10)), lwd=2, col="red", xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("EZ Deactivation (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

#plot(mtrates$denz_dom, mtrates$depth, type = "l", ylim=c(-3000, 0),  xlim=c(0,(max(mtrates$denz_dom)+max(mtrates$denz_dom)*0.10)), lwd=2, col="red")
#
#plottitle<-expression(paste("Mass Transfer Rates (mg ", m^-3," ", h^-1, ")", sep=""))
#title(plottitle, cex.main=1.5, line=0.2, outer=TRUE)

#dev.off()

}