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

linelist<-as.data.frame(seq(0, 0.05, 0.01))
#linelist<-as.data.frame(seq(0, 0.05, 0.005))
colnames(linelist)<-"maxe"

#---------------------
#     Interior
#---------------------
folderstate<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interior/StateVariables/"
statenames<-as.data.frame(list.files(folderstate))
statenames$maxe<-as.numeric(substr(statenames[,1], 20, 25))
hold<-merge(statenames, linelist)
statenames<-as.character(hold[,2])

folderbrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interior/BacteriaRates/"
bratesnames<-as.data.frame(list.files(folderbrates))
bratesnames$maxe<-as.numeric(substr(bratesnames[,1], 19, 24))
hold<-merge(bratesnames, linelist)
bratesnames<-as.character(hold[,2])

foldermrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interior/MassTransferRates/"
mratesnames<-as.data.frame(list.files(foldermrates))
mratesnames$maxe<-as.numeric(substr(mratesnames[,1], 23, 28))
hold<-merge(mratesnames, linelist)
mratesnames<-as.character(hold[,2])

analysis.ran<-as.data.frame(matrix(NA, length(statenames), 7))
colnames(analysis.ran)<-c("maxe", "qfrac.init", "qfrac.final", "ba.init", "ba.prod", "enz.prod", "h.prod")

for(a in 1:length(statenames)){

state<-read.table(paste(folderstate, statenames[a], sep=""))
colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

state<-state[362:378,]
brates<-brates[362:378,]
mrates<-mrates[362:378,]

analysis.ran$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
analysis.ran$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
analysis.ran$qfrac.final[a]<-state$ba[nrow(state)]/(quorum*pnum)
analysis.ran$ba.init[a]<-state$ba[1]
analysis.ran$ba.prod[a]<-sum(state$ba*brates$ba_grow)
analysis.ran$enz.prod[a]<-sum(state$ba*brates$epsilon)
analysis.ran$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

}

analysis.ran$qfrac.init<-analysis.ran$qfrac.init*100
analysis.ran$qfrac.final<-analysis.ran$qfrac.final*100

analysis.ran$ba.prod<-analysis.ran$ba.prod*1000  #convert mg to ug 
analysis.ran<-analysis.ran[order(analysis.ran$qfrac.init),]

all.interior<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonAll_Interior_LA_362to378min.txt", header=TRUE)
max.interior.baprod<-aggregate(all.interior$ba.prod, list(all.interior$qfrac.init), max)
colnames(max.interior.baprod)<-c("qfrac.init", "ba.prod")
max.interior.baprod$ba.prod<-max.interior.baprod$ba.prod/pnum*1000*1000  #convert to particle^-1 and from mg to ng 

#----------------------
#     Interception1
#----------------------
#folderstate<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception1/StateVariables/"
#statenames<-as.data.frame(list.files(folderstate))
#statenames$maxe<-as.numeric(substr(statenames[,1], 20, 25))
#hold<-merge(statenames, linelist)
#statenames<-as.character(hold[,2])

#folderbrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception1/BacteriaRates/"
#bratesnames<-as.data.frame(list.files(folderbrates))
#bratesnames$maxe<-as.numeric(substr(bratesnames[,1], 19, 24))
#hold<-merge(bratesnames, linelist)
#bratesnames<-as.character(hold[,2])

#foldermrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception1/MassTransferRates/"
#mratesnames<-as.data.frame(list.files(foldermrates))
#mratesnames$maxe<-as.numeric(substr(mratesnames[,1], 23, 28))
#hold<-merge(mratesnames, linelist)
#mratesnames<-as.character(hold[,2])

#analysis.int1<-as.data.frame(matrix(NA, length(statenames), 7))
#colnames(analysis.int1)<-c("maxe", "qfrac.init", "qfrac.final", "ba.init", "ba.prod", "enz.prod", "h.prod")

#for(a in 1:length(statenames)){

#state<-read.table(paste(folderstate, statenames[a], sep=""))
#colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

#brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
#colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

#mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
#colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

#state<-state[362:378,]
#brates<-brates[362:378,]
#mrates<-mrates[362:378,]

#analysis.int1$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
#analysis.int1$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
#analysis.int1$qfrac.final[a]<-state$ba[nrow(state)]/(quorum*pnum)
#analysis.int1$ba.init[a]<-state$ba[1]
#analysis.int1$ba.prod[a]<-sum(state$ba*brates$ba_grow)
#analysis.int1$enz.prod[a]<-sum(state$ba*brates$epsilon)
#analysis.int1$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

#}

#analysis.int1$qfrac.init<-analysis.int1$qfrac.init*100
#analysis.int1$qfrac.final<-analysis.int1$qfrac.final*100

#analysis.int1$ba.prod<-analysis.int1$ba.prod*1000  #convert mg to ug 
#analysis.int1<-analysis.int1[order(analysis.int1$qfrac.init),]


#----------------------
#     Interception2
#----------------------
folderstate<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception/StateVariables/"
statenames<-as.data.frame(list.files(folderstate))
statenames$maxe<-as.numeric(substr(statenames[,1], 20, 25))
hold<-merge(statenames, linelist)
statenames<-as.character(hold[,2])

folderbrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception/BacteriaRates/"
bratesnames<-as.data.frame(list.files(folderbrates))
bratesnames$maxe<-as.numeric(substr(bratesnames[,1], 19, 24))
hold<-merge(bratesnames, linelist)
bratesnames<-as.character(hold[,2])

foldermrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Interception/MassTransferRates/"
mratesnames<-as.data.frame(list.files(foldermrates))
mratesnames$maxe<-as.numeric(substr(mratesnames[,1], 23, 28))
hold<-merge(mratesnames, linelist)
mratesnames<-as.character(hold[,2])

analysis.int2<-as.data.frame(matrix(NA, length(statenames), 7))
colnames(analysis.int2)<-c("maxe", "qfrac.init", "qfrac.final", "ba.init", "ba.prod", "enz.prod", "h.prod")

for(a in 1:length(statenames)){

state<-read.table(paste(folderstate, statenames[a], sep=""))
colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

state<-state[362:378,]
brates<-brates[362:378,]
mrates<-mrates[362:378,]

analysis.int2$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
analysis.int2$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
analysis.int2$qfrac.final[a]<-state$ba[nrow(state)]/(quorum*pnum)
analysis.int2$ba.init[a]<-state$ba[1]
analysis.int2$ba.prod[a]<-sum(state$ba*brates$ba_grow)
analysis.int2$enz.prod[a]<-sum(state$ba*brates$epsilon)
analysis.int2$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

}

analysis.int2$qfrac.init<-analysis.int2$qfrac.init*100
analysis.int2$qfrac.final<-analysis.int2$qfrac.final*100

analysis.int2$ba.prod<-analysis.int2$ba.prod*1000  #convert mg to ug 
analysis.int2<-analysis.int2[order(analysis.int2$qfrac.init),]

all.interception<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonAll_Interception_LA_362to378min.txt", header=TRUE)
max.interception.baprod<-aggregate(all.interception$ba.prod, list(all.interception$qfrac.init), max)
colnames(max.interception.baprod)<-c("qfrac.init", "ba.prod")
max.interception.baprod$ba.prod<-max.interception.baprod$ba.prod/pnum*1000*1000  #convert to particle^-1 and from mg to ng 

#----------------------
#     Retention
#----------------------
folderstate<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Retention/StateVariables/"
statenames<-as.data.frame(list.files(folderstate))
statenames$maxe<-as.numeric(substr(statenames[,1], 20, 25))
hold<-merge(statenames, linelist)
statenames<-as.character(hold[,2])

folderbrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Retention/BacteriaRates/"
bratesnames<-as.data.frame(list.files(folderbrates))
bratesnames$maxe<-as.numeric(substr(bratesnames[,1], 19, 24))
hold<-merge(bratesnames, linelist)
bratesnames<-as.character(hold[,2])

foldermrates<-"/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/Retention/MassTransferRates/"
mratesnames<-as.data.frame(list.files(foldermrates))
mratesnames$maxe<-as.numeric(substr(mratesnames[,1], 23, 28))
hold<-merge(mratesnames, linelist)
mratesnames<-as.character(hold[,2])

analysis.ret<-as.data.frame(matrix(NA, length(statenames), 7))
colnames(analysis.ret)<-c("maxe", "qfrac.init", "qfrac.final", "ba.init", "ba.prod", "enz.prod", "h.prod")

for(a in 1:length(statenames)){

state<-read.table(paste(folderstate, statenames[a], sep=""))
colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

state<-state[362:378,]
brates<-brates[362:378,]
mrates<-mrates[362:378,]

analysis.ret$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
analysis.ret$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
analysis.ret$qfrac.final[a]<-state$ba[nrow(state)]/(quorum*pnum)
analysis.ret$ba.init[a]<-state$ba[1]
analysis.ret$ba.prod[a]<-sum(state$ba*brates$ba_grow)
analysis.ret$enz.prod[a]<-sum(state$ba*brates$epsilon)
analysis.ret$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

}

analysis.ret$qfrac.init<-analysis.ret$qfrac.init*100
analysis.ret$qfrac.final<-analysis.ret$qfrac.final*100

analysis.ret$ba.prod<-analysis.ret$ba.prod*1000  #convert mg to ug 
analysis.ret<-analysis.ret[order(analysis.ret$qfrac.init),]

all.retention<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonAll_Retention_LA_362to378min.txt", header=TRUE)
max.retention.baprod<-aggregate(all.retention$ba.prod, list(all.retention$qfrac.init), max)
colnames(max.retention.baprod)<-c("qfrac.init", "ba.prod")
max.retention.baprod$ba.prod<-max.retention.baprod$ba.prod/pnum*1000*1000  #convert to particle^-1 and from mg to ng 

#------------------------------------
#     Figure 
#------------------------------------
library(fields)

OutGraph<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Graphs/Figure_Scenarios5.ps", sep="")
postscript(OutGraph, family="Helvetica", width=3, height=8, pointsize=15)

#quartz(width=3, height=8)
par(oma=c(7,3,1,1))
par(mar=c(1,2,0.3,0.5)+0.1)
par(mfrow=c(3,1))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)
#col="lightgrey"

# Interior Arrangement
listmaxe<-unique(analysis.ran$maxe)
toplot<-subset(analysis.ran, analysis.ran$maxe==listmaxe[1])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
plot(toplot$qfrac.init, toplot$ba.prod, type="l", ylim=c(-20, 60), xlim=c(0,110), col=tim.colors()[1], xlab="", ylab="", xaxt="n", las=1, lwd=2)
abline(h=0, lwd=1)
#mtext(expression(paste("Bacterial production  (", mu, g, " ", m^{-3}, " ", h^{-1}, ")")), side=2, line=2.5, cex=0.8)
axis(side=1, labels=FALSE)

for(i in 2:length(listmaxe)){
toplot<-subset(analysis.ran, analysis.ran$maxe==listmaxe[i])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
lines(toplot$qfrac.init, toplot$ba.prod, col=tim.colors()[i*10], lwd=2)
	}

#lines(max.interior.baprod$qfrac.init, max.interior.baprod$ba.prod, lty=3, lwd=2, col="black")

legend(58, 0.98*60, sprintf("%1.2f", listmaxe), seg.len=0.60, ncol=2, lwd=2, cex=1, xjust=0, col=tim.colors()[c(1,20,30,40,50,60)], y.intersp=0.9, x.intersp=0.8, title=expression(epsilon))



#Interception 
listmaxe<-unique(analysis.int2$maxe)
toplot<-subset(analysis.int2, analysis.int2$maxe==listmaxe[1])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
plot(toplot$qfrac.init, toplot$ba.prod, type="l", ylim=c(-20, 60), xlim=c(0,110), col=tim.colors()[1], xlab="", ylab="", xaxt="n", yaxt="n",las=1, lwd=2)
abline(h=0, lwd=1)
#mtext(expression(paste("Bacterial production  (", mu, g, " ", m^{-3}, " ", h^{-1}, ")")), side=2, line=2.5, cex=0.8)
axis(side=1, labels=FALSE)
axis(side=2, labels=TRUE, las=2)

for(i in 2:length(listmaxe)){
toplot<-subset(analysis.int2, analysis.int2$maxe==listmaxe[i])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
lines(toplot$qfrac.init, toplot$ba.prod, col=tim.colors()[i*10], lwd=2)
	}	
#lines(max.interception.baprod$qfrac.init, max.interception.baprod$ba.prod, lty=3, lwd=2, col="black")	
	
#Retention 
listmaxe<-unique(analysis.ret$maxe)
toplot<-subset(analysis.ret, analysis.ret$maxe==listmaxe[1])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
plot(toplot$qfrac.init, toplot$ba.prod, type="l", ylim=c(-20, 200), xlim=c(0,110), col=tim.colors()[1], xlab="", ylab="", xaxt="n", yaxt="n",las=1, lwd=2)
abline(h=0, lwd=1)
axis(side=1, labels=TRUE)
axis(side=2, labels=TRUE, las=2)

per<-seq(0, 100, 20)
babund<-round(((per/100)*quorum/cperb)/1E6, digits=0)
axis(side=1, col="black", line=4, labels=babund, at=per)

for(i in 2:length(listmaxe)){
toplot<-subset(analysis.ret, analysis.ret$maxe==listmaxe[i])
toplot$ba.prod<-toplot$ba.prod/pnum*1000 #convert to particle^-1 and from ug to ng
lines(toplot$qfrac.init, toplot$ba.prod, col=tim.colors()[i*10], lwd=2)
	}		
#lines(max.retention.baprod$qfrac.init, max.retention.baprod$ba.prod, lty=3, lwd=2, col="black")	
	
mtext("Bacterial coverage of particle (%)", side=1, line=2.4, cex=0.8, adj=0.5)
mtext(expression(paste("Bacterial abundance (",10^{6}, " ", particle^{-1}, ")")), side=1, line=6.4, cex=0.8, adj=0.5)
mtext(expression(paste("Bacterial production (", ng, " C ", particle^{-1}, " ", h^{-1}, ")")), side=2, line=0.7, cex=0.8, outer=TRUE)

dev.off()
