Folder<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/CFlux1D/Interior/StateVariables/")
Folder2<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/CFlux1D/Interception/StateVariables/")
Folder3<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/CFlux1D/Retention/StateVariables/")

statefiles<-list.files(Folder)
statefiles2<-list.files(Folder2)
statefiles3<-list.files(Folder3)

#statefiles<-subset(statefiles, statefiles=="StateVariables_Sink2.00_Size0.0010_Qfrac1.00.out")
#statefiles2<-subset(statefiles2, statefiles2=="StateVariables_Sink2.00_Size0.0010_Qfrac1.00.out")
#statefiles3<-subset(statefiles3, statefiles3=="StateVariables_Sink2.00_Size0.0010_Qfrac1.00.out")

#----------------------
#	INTERIOR CFLUX
#----------------------

cflux.interior<-as.data.frame(matrix(NA, 386, (1+2*length(statefiles))))
colnames(cflux.interior)[1]<-"depth"
cflux.interior$depth<-seq(150,4000,10)

count<-2
for(a in 1:length(statefiles)){

state<-read.table(paste(Folder, statefiles[a], sep=""))
colnames(state)<-c("depth", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom", "pomflux", "baflux", "enzflux", "hflux", "denzflux")
state$flux<-state$pomflux+state$baflux
state$flux<-state$flux*24
state$flux2<-state$pomflux+state$baflux+state$enzflux+state$hflux+state$denzflux
state$flux2<-state$flux2*24
#state3<-merge(state3, depth2, all.y=TRUE)
state$depth<-abs(state$depth)*-1
#print(sum((state$flux-min(state$flux))*state$depth)/sum(state$flux-min(state$flux)))
#print(sum((state$flux2-min(state$flux2))*state$depth)/sum(state$flux2-min(state$flux2)))


fluxlist<-as.data.frame(seq(2, max(state$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state[,c("depth", "flux")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepth<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepth[1])

fluxlist<-as.data.frame(seq(2, max(state$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state[,c("depth", "flux2")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepthb<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepthb[1])

colnames(cflux.interior)[count:(count+1)]<-c(paste(substr(statefiles[a],36,44),"_Cp",sep=""),paste(substr(statefiles[a],36,44),"_Cs",sep=""))

cflux.interior[,count]<-state$flux
cflux.interior[,count+1]<-state$flux2
count<-count+2

}

#----------------------
#	INTERCEPTION CFLUX
#----------------------

cflux.interception<-as.data.frame(matrix(NA, 386, (1+2*length(statefiles2))))
colnames(cflux.interception)[1]<-"depth"
cflux.interception$depth<-seq(150,4000,10)

count<-2
for(b in 1:length(statefiles2)){

state2<-read.table(paste(Folder2, statefiles2[b], sep=""))
colnames(state2)<-c("depth", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom", "pomflux", "baflux", "enzflux", "hflux", "denzflux")
state2$flux<-state2$pomflux+state2$baflux
state2$flux<-state2$flux*24
state2$flux2<-state2$pomflux+state2$baflux+state2$enzflux+state2$hflux+state2$denzflux
state2$flux2<-state2$flux2*24
#state3<-merge(state3, depth2, all.y=TRUE)
state2$depth<-abs(state2$depth)*-1
#print(sum((state2$flux-min(state2$flux))*state2$depth)/sum(state2$flux-min(state2$flux)))
#print(sum((state2$flux2-min(state2$flux2))*state2$depth)/sum(state2$flux2-min(state2$flux2)))

fluxlist<-as.data.frame(seq(2, max(state2$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state2[,c("depth", "flux")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepth2<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepth2[1])

fluxlist<-as.data.frame(seq(2, max(state2$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state2[,c("depth", "flux2")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepth2b<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepth2b[1])

colnames(cflux.interception)[count:(count+1)]<-c(paste(substr(statefiles2[b],36,44),"_Cp",sep=""),paste(substr(statefiles2[b],36,44),"_Cs",sep=""))

cflux.interception[,count]<-state2$flux
cflux.interception[,count+1]<-state2$flux2
count<-count+2

}

#----------------------
#	RETENTION CFLUX
#----------------------

cflux.retention<-as.data.frame(matrix(NA, 386, (1+2*length(statefiles))))
colnames(cflux.retention)[1]<-"depth"
cflux.retention$depth<-seq(150,4000,10)

count<-2
for(z in 1:length(statefiles3)){

state3<-read.table(paste(Folder3, statefiles3[z], sep=""))
colnames(state3)<-c("depth", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom", "pomflux", "baflux", "enzflux", "hflux", "denzflux")
state3$flux<-state3$pomflux+state3$baflux
state3$flux<-state3$flux*24
state3$flux2<-state3$pomflux+state3$baflux+state3$enzflux+state3$hflux+state3$denzflux
state3$flux2<-state3$flux2*24
#state3<-merge(state3, depth2, all.y=TRUE)
state3$depth<-abs(state3$depth)*-1
#print(sum((state3$flux-min(state3$flux))*state3$depth)/sum(state3$flux-min(state3$flux)))
#print(sum((state3$flux2-min(state3$flux2))*state3$depth)/sum(state3$flux2-min(state3$flux2)))

fluxlist<-as.data.frame(seq(2, max(state3$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state3[,c("depth", "flux")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepth3<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepth3[1])

fluxlist<-as.data.frame(seq(2, max(state3$flux), 0.1))
colnames(fluxlist)<-"flux"
fluxtable<-state3[,c("depth", "flux2")]
fluxtable$flux<-round(fluxtable$flux, 1)
fluxtable$flux<-as.numeric(fluxtable$flux)
fluxtable$depth<-abs(fluxtable$depth)
fluxfind<-merge(fluxtable, fluxlist, all=TRUE)
x<-fluxfind$flux
y<-fluxfind$depth
fd<-approx(x,y,x, rule=1)
fluxfind$depth<-fd$y
fluxfind<-fluxfind[order(fluxfind$depth),]
fluxfind<-unique(fluxfind)
rownames(fluxfind)<-c()
fluxfind$flux<-round(fluxfind$flux, digits=1)
efoldflux<-round((1/exp(1))*fluxfind$flux[1], digits=1)
efolddepth3b<-round(fluxfind$depth[which(fluxfind$flux==efoldflux)], digits=0)
print(efolddepth3b[1])

colnames(cflux.retention)[count:(count+1)]<-c(paste(substr(statefiles3[z],36,44),"_Cp",sep=""),paste(substr(statefiles3[z],36,44),"_Cs",sep=""))

cflux.retention[,count]<-state3$flux
cflux.retention[,count+1]<-state3$flux2
count<-count+2

}

#----------------------------------------------------
#  Calculate Mean Remineralization E-folding Depth
#----------------------------------------------------

#edepthlist<-do.call("c", lapply(statefiles, as.name))

#print(round(mean(edepthlist), digits=0))
#print(round(sd(edepthlist), digits=0))

#---------------------------
#  Analyze Trap Data
#---------------------------

Trapdata<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/BATS_trapC.txt", header=TRUE)
Trapdata<-subset(Trapdata, Trapdata$Cavg!=-9.99)
Trapdata<-subset(Trapdata, Trapdata$yymmdd1>=19900000)
Trapdata$Depth<-abs(Trapdata$Depth)*-1
Trapdata<-subset(Trapdata, Trapdata$Depth!=-400)

trap150<-subset(Trapdata, Trapdata$Depth==-150)
mean150<-round(mean(trap150$Cavg, na.rm=TRUE), digits=1)
sd150<-sd(trap150$Cavg, na.rm=TRUE)
nearmean<-subset(trap150, trap150$Cavg>26.8 & trap150$Cavg<28.8)

nearmean<-nearmean[,c("X.Cruise", "yymmdd1", "yymmdd2", "Lat1", "Lat2", "Long1", "Long2")]
Trapdata<-merge(Trapdata, nearmean)

trap200<-subset(Trapdata, Trapdata$Depth==-200)
mean200<-mean(trap200$Cavg, na.rm=TRUE)
sd200<-sd(trap200$Cavg, na.rm=TRUE)

trap300<-subset(Trapdata, Trapdata$Depth==-300)
mean300<-mean(trap300$Cavg, na.rm=TRUE)
sd300<-sd(trap300$Cavg, na.rm=TRUE)

#---------------------------
#   Calculate Martin Curve
#---------------------------
cavg<-aggregate(Trapdata, list(Trapdata$Depth), mean)

depth<-c(150, 200, 300, 500, 1500, 3200)
cflux<-c(cavg$Cavg[which(cavg$Depth==-150)], cavg$Cavg[which(cavg$Depth==-200)],cavg$Cavg[which(cavg$Depth==-300)], 4.1, 2.4, 1.7)

#plot(log(depth), log(cflux))

#summary(lm(log(cflux)~log(depth)))

mdepth<-c(seq(150, 600, 50), seq(700, 4000, 100))*-1
martin<-cflux[1]*(abs(state3$depth)/150)^(-0.94)




#------------------------------------------------------
#   Calculate Martin Remineralization E-folding Depth
#------------------------------------------------------

mfluxlist<-as.data.frame(seq(2, cflux[1], 0.1))
colnames(mfluxlist)<-"mflux"
mfluxtable<-as.data.frame(cbind(state3$depth, martin))
colnames(mfluxtable)<-c("depth", "mflux")
mfluxtable$mflux<-round(mfluxtable$mflux, 1)
mfluxtable$mflux<-as.numeric(mfluxtable$mflux)
mfluxtable$depth<-abs(mfluxtable$depth)
mfluxfind<-merge(mfluxtable, mfluxlist, all=TRUE)
x<-mfluxfind$mflux
y<-mfluxfind$depth
mfd<-approx(x,y,x, rule=2)
mfluxfind$depth<-mfd$y
mfluxfind<-mfluxfind[order(mfluxfind$depth),]
rownames(mfluxfind)<-c()
mfluxfind<-unique(mfluxfind)
mefoldflux<-round((1/exp(1))*mfluxfind$mflux[1], digits=1)
mefolddepth<-round(mfluxfind$depth[which(mfluxfind$mflux==mefoldflux)], digits=0)




#------------------------------
#   MAKE PLOT
#------------------------------
OutGraph<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Graphs/Figure_CFlux2.ps", sep="")
postscript(OutGraph, family="Helvetica", width=6, height=3, pointsize=12)

Qfraclist<-c("Qfrac0.25_Cp", "Qfrac0.50_Cp", "Qfrac1.00_Cp")

depth<-as.data.frame(c(seq(150, 600, 50), seq(700, 4000, 100)))
colnames(depth)<-"depth"

cflux.interior<-merge(cflux.interior, depth, all.y=TRUE)
cflux.interception<-merge(cflux.interception, depth, all.y=TRUE)

#quartz(width=6,height=3)
par(mfrow=c(1,3))
par(oma=c(1,3,3,1))
par(mar=c(0.5,1.5,1.5,0)+0.1)
par(xaxs="i")
par(yaxs="i")

library(RColorBrewer)
colpick<-brewer.pal(4, "Spectral")
#colpick<-brewer.pal(4, "RdYlBu")
#colpick<-brewer.pal(4, "Dark2")

for(j in 1:length(Qfraclist)){

p<-which(colnames(cflux.interior)==Qfraclist[j])
q<-which(colnames(cflux.interception)==Qfraclist[j])
r<-which(colnames(cflux.retention)==Qfraclist[j])
plot(cflux.interior[,p], abs(cflux.interior$depth)*-1, ylim=c(-4000,0), ylab="", xlab="", xlim=c(0,30), xaxt="n", yaxt="n", col=colpick[1], type="p", pch=17, cex=1)
points(cflux.interception[,q], abs(cflux.interception$depth)*-1, col=colpick[2], pch=16, cex=0.80)
lines(cflux.retention[,r], abs(cflux.retention$depth)*-1, lwd=3.5, col=colpick[4], lty=1)
lines(cflux.retention[,r+1], abs(cflux.retention$depth)*-1, lwd=2, col="cyan2", lty=1)
axis(side=3, tick=TRUE, cex.axis=1.2, labels=c(0,"", 10,"", 20, "", 30), at=c(0, 5, 10, 15, 20, 25, 30))
if(j==1){axis(side=2, las=2, tick=TRUE, labels=c(0,1000,2000,3000,4000), at=c(0,-1000,-2000,-3000,-4000), cex.axis=1.2)
mtext("Depth (m)", side=2, line=3, cex=1)}else{axis(side=2, las=2, tick=TRUE, labels=FALSE, at=c(0,-1000,-2000,-3000,-4000))}

lines(martin, abs(state3$depth)*-1, col="grey60", lty=2, lwd=1.5)

#arrows(mean150, -150, (mean150+sd150), -150, length=0.05, angle=90, lwd=1, col="black")
#arrows(mean150, -150, (mean150-sd150), -150, length=0.05, angle=90, lwd=1, col="black")
arrows(mean200, -200, (mean200+sd200), -200, length=0.05, angle=90, lwd=1, col="black")
arrows(mean200, -200, (mean200-sd200), -200, length=0.05, angle=90, lwd=1, col="black")
arrows(mean300, -300, (mean300+sd300), -300, length=0.05, angle=90, lwd=1, col="black")
arrows(mean300, -300, (mean300-sd300), -300, length=0.05, angle=90, lwd=1, col="black")
arrows(4.1, -500, (4.1+3.1), -500, length=0.05, angle=90, lwd=1, col="black")
arrows(4.1, -500, (4.1-3.1), -500, length=0.05, angle=90, lwd=1, col="black")
arrows(2.4, -1500, (2.4+1.3), -1500, length=0.05, angle=90, lwd=1, col="black")
arrows(2.4, -1500, (2.4-1.3), -1500, length=0.05, angle=90, lwd=1, col="black")
arrows(1.7, -3200, (1.7+0.78), -3200, length=0.05, angle=90, lwd=1, col="black")
arrows(1.7, -3200, (1.7-0.78), -3200, length=0.05, angle=90, lwd=1, col="black")

#points(mean150, -150, pch=16, col="black", cex=0.6)
points(mean200, -200, pch=18, col="black", cex=1.6)
points(mean300, -300, pch=18, col="black", cex=1.6)
points(4.1, -500, pch=18, col="black", cex=1.6)
points(2.4, -1500, pch=18, col="black", cex=1.6)
points(1.7, -3200, pch=18, col="black", cex=1.6)


}

legend(7, -2200, box.lwd=0, seg.len=1, cex=1.2, ncol=1, c( "Interior","Interception","Retention", "Retention+","Sediment Trap", "Martin Curve"), lwd=c(3.5, 3.5, 3.5, 2, 3.5, 3.5), col=c(colpick[1], colpick[2], colpick[4], "cyan", "black", "grey60"), y.intersp=0.9)

xaxis<-expression(paste("Carbon flux (mg ", m^-2," ", d^-1, ")", sep=""))
mtext(xaxis, side=3, line=1, cex=1, outer=TRUE)

dev.off()


