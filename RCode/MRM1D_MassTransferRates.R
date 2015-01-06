#Change working directory to MicrobeReminModel_v1.0/RCode using the setwd() command or using the Misc menu.  
#getwd() #shows the current working directory for the R gui.

#----------------------
# Read in Model Output
#----------------------

mtfolder<-paste("../Output/MRM1D/MassTransferRates/")
mtfiles<-list.files(mtfolder)

for(i in 1:length(mtfiles)){

mtrates<-read.table(paste(mtfolder, mtfiles[i], sep=""))
colnames(mtrates)<-c("depth", "pnum", "pom_dgd", "h_dom", "ehl_k", "enz_dom", "denz_dom", "bacover") 
mtrates$depth<-mtrates$depth*-1

#----------------
#  Create Plot
#----------------

outname<-substr(mtfiles[i], 1, nchar(mtfiles[i])-4)
OutGraph<-paste("../Graphs/MRM1D1_", outname, ".ps", sep="")
postscript(OutGraph, family="Helvetica", width=9, height=3)

#quartz(width=7.5, height=6)
par(oma=c(2,4,2,2))
par(mar=c(1,1,4,1))
par(mfrow=c(1,5))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)

pdepth<--4000

#-------------------------
# Plot Mass Transfer Rates
#--------------------------

plot(mtrates$pom_dgd, mtrates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(mtrates$pom_dgd)+max(mtrates$pom_dgd)*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="Depth (m)", xlab="")
title(expression(paste("Hydrolysis of POM (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, las=2, tick=TRUE, labels=seq(0, -1*pdepth, 1000), at=seq(0,pdepth,-1000))

plot(mtrates$ehl_k, mtrates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(mtrates$ehl_k)+max(mtrates$ehl_k)*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("EP Deactivation (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(mtrates$h_dom, mtrates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(mtrates$h_dom)+max(mtrates$h_dom)*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Mass Transfer HP (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(mtrates$enz_dom, mtrates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(mtrates$enz_dom)+max(mtrates$enz_dom)*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Mass Transfer EP (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)

plot(mtrates$denz_dom, mtrates$depth, type = "l", ylim=c(pdepth, 0),  xlim=c(0,(max(mtrates$denz_dom)+max(mtrates$denz_dom)*0.10)), lwd=2, xaxt="n", yaxt="n", ylab="", xlab="")
title(expression(paste("Mass Transfer XP (", h^-1, ")", sep="")), line=2.7)
axis(side=3)
axis(side=2, labels=FALSE)


#plottitle<-expression(paste("Mass Transfer Rates (mg ", m^-3," ", h^-1, ")", sep=""))
#title(plottitle, cex.main=1.5, line=0.2, outer=TRUE)

dev.off()

}