##Data Input
library(caTools)
tempdata<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/CTDTables_PP/BATS_CTDTemp_10m.txt", header=TRUE)
tempdata<-tempdata[,c("Depth", "mean")]
#plot(tempdata$mean, tempdata$Depth*-1, type="l")


Sdata<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/CTDTables_PP/BATS_CTDSal_10m.txt", header=TRUE)
Sdata<-Sdata[,c("Depth", "mean")]
#plot(Sdata$mean, Sdata$Depth*-1, type="l")


Rhodata<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/CTDTables_PP/BATS_CTDRho_10m.txt", header=TRUE)
Rhodata<-Rhodata[,c("Depth", "mean")]
Rhodata$mean<-Rhodata$mean+1000
#plot(Rhodata$mean, Rhodata$Depth*-1, type="l")


DOCdata<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/BATS_DOC1990-2010.txt", header=TRUE)
DOCdata<-subset(DOCdata, DOCdata$DOC>=0)
DOCdata$DOC<-DOCdata$DOC*1.02632 #convert from umol/kg to mmol/m^3
DOCdata$DOC<-(DOCdata$DOC-43.6)*12 #remove refractory component, convert from mmol/m^3 to mg/m^3
DOCdata$DOC<-ifelse(DOCdata$DOC<0, 0, DOCdata$DOC)

depth<-seq(0, 4000, 10)

## Semi-labile DOC profile ##
OceanDOC<-DOCdata[,c("Depth", "DOC")]
docp<-as.data.frame(matrix(NA, length(depth), 2))
colnames(docp)<-c("Depth", "doc")
docp[,1]<-depth
for(z in 1:nrow(docp)){
	upper<-docp$Depth[z]-5
	lower<-docp$Depth[z]+5
	hold<-subset(OceanDOC, OceanDOC$Depth>=upper & OceanDOC$Depth<=lower)
	if(nrow(hold)>100){
	print(nrow(hold))
	docp$doc[z]<-mean(hold$DOC, na.rm=TRUE)}
	}

docp2<-subset(docp, is.na(docp$doc)==FALSE)
docp2<-subset(docp2, docp2$Depth<1500)
lo<-loess(docp2$doc~docp2$Depth)
docp2$lo<-predict(lo)

docp<-merge(docp, docp2, all=TRUE)

xy<-approx(docp$Depth, docp$lo, docp$Depth, rule=2)
docp$lo<-xy$y
#docp<-subset(docp, docp$Depth>=150)



xl <- seq(min(docp$Depth),max(docp$Depth), (max(docp$Depth) - min(docp$Depth))/10000)

smore<-predict(lo,xl)

hold<-as.data.frame(cbind(xl, smore))
colnames(hold)<-c("Depth", "doc2")
docp3<-merge(docp, hold)

diffhold<-diff(docp3$doc2)
docp3$diff<-rep(NA)
docp3$diff[2:nrow(docp3)]<-diffhold

docp3$doc2<-ifelse(docp3$diff>0, NA, docp3$doc2)

docp<-subset(docp3, docp3$Depth>=150)

mindoc2<-min(docp$doc2, na.rm=TRUE)

docp$doc2<-ifelse(is.na(docp$doc2)==TRUE, mindoc2, docp$doc2)


docp$doc2<-round(docp$doc2, digits=2)

plot(docp$doc2, docp$Depth*-1, type="l")



##Final Synthesis##

makeprofile<-cbind(tempdata, Sdata$mean, Rhodata$mean, docp$doc2)
colnames(makeprofile)<-c("Depth", "temp", "sal", "rho", "sldoc")

#plot(makeprofile$sldoc, makeprofile$Depth*-1, type="l")

outfile<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/TrapData/BATSin/Profiles/Profile_Mean10m.in", sep="")

write.table(makeprofile, outfile, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")


#Black and white version
#OutGraph<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/TrapData/BATSgraphs/InputProfilesBlackWhite.ps", sep="")
#postscript(OutGraph, family="Helvetica", width=4, height=8)
quartz(width=4, height=7)
#par(mfrow=c(1,2))
#from http://www.r-bloggers.com/multiple-y-axis-in-a-r-plot/
par(oma=c(2,2,2,2))
par(mar=c(1, 4 , 14, 1) + 0.1)
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=0.8)
#par(mar=c(3, 5 , 14, 3) + 0.1)
plot(makeprofile$temp, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="Depth (m)",type="l",col="black", main="",xlim=c(0,20), lwd=2)
axis(2, col="black",lwd=1, las=2, tick=TRUE, labels=c(0,1000,2000,3000,4000), at=c(0,-1000,-2000,-3000,-4000))
axis(3, col="black",lwd=1)
mtext(3,text="Temperature (°C)",line=2)
#mtext(2,text="Depth (m)", line=2.5)

par(new=T)
plot(makeprofile$sal, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="black", main="",xlim=c(34,37), lwd=2, lty=2)
axis(3, col="black",lwd=1, line=3.5)
mtext(3,text="Salinity (psu)",line=5.5)

par(new=T)
plot(makeprofile$rho, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="black", main="",xlim=c(1026,1028), lwd=2, lty=3)
axis(3, col="black",lwd=1, line=7)
rhotext<-expression(paste("Density (kg ", m^-3, ")", sep=""))
mtext(3,text=rhotext,line=9)

#quartz(width=4, height=7)
par(new=T)
plot(makeprofile$sldoc, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="grey50", main="",xlim=c(0,200), lwd=2)
#axis(2, col="black",lwd=2)
axis(3, col="black",lwd=1, line=10.5)
sltext<-expression(paste("Semi-labile DOM (mg ", m^-3, ")", sep=""))
mtext(3,text=sltext,line=12.5)
#mtext(2,text="Depth (m)", line=2.5)

#dev.off()



##Color version
#OutGraph<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/TrapData/BATSgraphs/InputProfilesColor.ps", sep="")
#postscript(OutGraph, family="Helvetica", width=4, height=8)
##quartz(width=4, height=8)
##par(mfrow=c(1,2))
##from http://www.r-bloggers.com/multiple-y-axis-in-a-r-plot/
#par(mar=c(7, 5 , 14, 3) + 0.1)
#plot(makeprofile$temp, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="red", main="",xlim=c(0,20), lwd=2)
#axis(2, col="black",lwd=2)
#axis(3, col="red",lwd=2)
#mtext(3,text="Temperature (°C)",line=2)
#mtext(2,text="Depth (m)", line=2.5)
#
#par(new=T)
#plot(makeprofile$sal, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="green", main="",xlim=c(34,37), lwd=2)
#axis(3, col="green",lwd=2, line=3.5)
#mtext(3,text="Salinity (psu)",line=5.5)
#
#par(new=T)
#plot(makeprofile$rho, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="blue", main="",xlim=c(1026,1028), lwd=2)
#axis(3, col="blue",lwd=2, line=7)
#mtext(3,text="Density (kg m^-3)",line=9)
#
##quartz(width=4, height=7)
#par(new=T)
#plot(makeprofile$sldoc, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="grey50", main="",xlim=c(0,200), lwd=2)
##axis(2, col="black",lwd=2)
#axis(3, col="grey50",lwd=2, line=10.5)
#mtext(3,text="Semi-labile DOM (mg m^-3)",line=12.5)
##mtext(2,text="Depth (m)", line=2.5)
#
#dev.off()
#
##Color version2
##OutGraph<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/TrapData/BATSgraphs/InputProfilesColor.ps", sep="")
##postscript(OutGraph, family="Helvetica", width=4, height=8)
#quartz(width=4, height=8)
##par(mfrow=c(1,2))
##from http://www.r-bloggers.com/multiple-y-axis-in-a-r-plot/
##par(mar=c(3, 5 , 14, 3) + 0.1)
#par(mar=c(7, 5 , 7, 3) + 0.1)
#plot(makeprofile$temp, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="red", main="",xlim=c(0,20), lwd=2)
#axis(2, col="black",lwd=2)
#axis(3, col="red",lwd=2)
#mtext(3,text="Temperature (°C)",line=2)
#mtext(2,text="Depth (m)", line=2.5)
#
#par(new=T)
#plot(makeprofile$sal, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="green", main="",xlim=c(34,37), lwd=2)
#axis(3, col="green",lwd=2, line=3.5)
#mtext(3,text="Salinity (psu)",line=5.5)
#
#par(new=T)
#plot(makeprofile$rho, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="blue", main="",xlim=c(1026,1028), lwd=2)
#axis(1, col="blue",lwd=2)
#mtext(1,text="Density (kg m^-3)",line=2)
#
##quartz(width=4, height=7)
#par(new=T)
#plot(makeprofile$sldoc, makeprofile$Depth*-1, axes=F, ylim=c(-4000,0), xlab="", ylab="",type="l",col="grey50", main="",xlim=c(0,200), lwd=2)
##axis(2, col="black",lwd=2)
#axis(1, col="grey50",lwd=2, line=3.5)
#mtext(1,text="Semi-labile DOM (mg m^-3)",line=5.5)
##mtext(2,text="Depth (m)", line=2.5)
#
#dev.off()
#
#
