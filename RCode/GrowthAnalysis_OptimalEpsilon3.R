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
scenario = "Retention"

filepath = paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/Growth3/", scenario, "/", sep="")

#tsegs1<-seq(109,195,9)
#tsegs2<-seq(117,200,9)

#tsegs1<-seq(181,185,9)
#tsegs2<-seq(189,190,9)

tsegs1<-seq(362,365,18)
tsegs2<-seq(378,380,18)


#--------------------------
#     Find Epsilon Optimum
#--------------------------
folderstate<-paste(filepath, "StateVariables/", sep="")
statenames<-list.files(folderstate)

folderbrates<-paste(filepath, "BacteriaRates/", sep="")
bratesnames<-list.files(folderbrates)

foldermrates<-paste(filepath, "MassTransferRates/", sep="")
mratesnames<-list.files(foldermrates)


for(b in 1:length(tsegs1)){

analysis.none<-as.data.frame(matrix(NA, length(statenames), 7))
colnames(analysis.none)<-c("maxe", "qfrac.init", "qfrac.final", "ba.init", "ba.prod", "enz.prod", "h.prod")

for(a in 1:length(statenames)){

state<-read.table(paste(folderstate, statenames[a], sep=""))
colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

state<-state[tsegs1[b]:tsegs2[b],]
brates<-brates[tsegs1[b]:tsegs2[b],]
mrates<-mrates[tsegs1[b]:tsegs2[b],]

analysis.none$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
analysis.none$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
analysis.none$qfrac.final[a]<-state$ba[nrow(state)]/(quorum*pnum)
analysis.none$ba.init[a]<-state$ba[1]
analysis.none$ba.prod[a]<-sum(state$ba*brates$ba_grow)
analysis.none$enz.prod[a]<-sum(state$ba*brates$epsilon)
analysis.none$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

}

analysis.none$qfrac.init<-analysis.none$qfrac.init*100
analysis.none$qfrac.final<-analysis.none$qfrac.final*100

opt.none<-aggregate(analysis.none$ba.prod, list(analysis.none$qfrac.init), max)
colnames(opt.none)<-c("qfrac.init", "ba.prod")
opt.none<-merge(opt.none, analysis.none)
opt.none<-opt.none[order(opt.none$qfrac.init),]

tname1<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonAll_", scenario, "_LA_", sprintf("%03.f", tsegs1[b]), "to", sprintf("%03.f", tsegs2[b]), "min.txt", sep="")
write.table(analysis.none, tname1, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


#tname1<-paste("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_", scenario, "_LA_", sprintf("%03.f", tsegs1[b]), "to", sprintf("%03.f", tsegs2[b]), "min.txt", sep="")
#write.table(opt.none, tname1, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

}

#------------------------------------
#     Write Retention Lookup Table
#------------------------------------
opt.none<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Retention_LA_362to378min.txt", header=TRUE)
opt.none2<-opt.none[,c("qfrac.init", "maxe")]
opt.none2<-subset(opt.none2, opt.none2$qfrac.init <= 100)
tname2<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Fortran/RetentionEpsilon2.txt", sep="")
write.table(opt.none2, tname2, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")

#------------------------------------
#     Figure Optimal Epsilon Max
#------------------------------------

quartz(width=4, height=3.25)
par(oma=c(3,3,1,1))
par(mar=c(1.5,1,1,0.5)+0.1)
par(mfrow=c(1,1))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)
#col="lightgrey"
# Interior Arrangement

#folder<-"/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/Spinup_Calc_Retention/"
#filelist<-list.files(folder)

#opt.none<-read.table(paste(folder, filelist[1], sep=""), header=TRUE)

#opt.none$maxe[which(opt.none$qfrac.init==100)]<-opt.none$maxe[110]
#opt.none$maxe[which(opt.none$qfrac.init==101)]<-opt.none$maxe[110]

#opt.none<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Retention_100to110min.txt", header=TRUE)
plot(opt.none$qfrac.init, opt.none$maxe, type="l", xlim=c(0,110), ylim=c(0,0.05), xaxt="n", lwd=2)
axis(side=1, labels=TRUE)

#for(d in 2:length(filelist)){
#d<-10
#opt.none<-read.table(paste(folder, filelist[d], sep=""), header=TRUE)
#lines(opt.none$qfrac.init, opt.none$maxe)
#}

#lines(opt.none$qfrac.init, opt.none$maxe, col="blue")

mtext("Initial bacterial coverage of particle (%)", side=1, cex=1, line=0.7 , adj=0.5, outer=TRUE)
mtext(expression(paste("Optimal  ", epsilon)), side=2, cex=1, line=1.2 , adj=0.5, outer=TRUE)


#dev.off()


#-------------------------------------------------
#    Calculate Exponent for QS Equation
#-------------------------------------------------

#EpsilonOpt<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Interior_LA_362to378min.txt", header=TRUE)
#emax=0.027
#et = 37
#ke = 52


EpsilonOpt<-read.table("/Users/kasmith/Dropbox/Manuscript/QuorumSensingPaper/Results/EpsilonOpt_Interception_LA_362to378min.txt", header=TRUE)
emax=0.025
et = 31
ke = 43


x<-seq(1, 3, 0.01)

R2<-as.data.frame(matrix(NA, length(x), 2))
colnames(R2)<-c("xvalue", "R2value")

for(g in 1:length(x)){

qfrac<-seq(1, 110, 1)
ep<-emax*abs(qfrac-et)^x[g]/((ke-et)^x[g]+abs(qfrac-et)^x[g])
ep2<-c(diff(ep),1)
ep<-ifelse(ep2<0, 0, ep)

R2$xvalue[g]<-x[g]
R2$R2value[g]<-summary(lm(EpsilonOpt$maxe~ep))$adj.r.squared


}

R2[which(R2$R2value==max(R2$R2value)),]