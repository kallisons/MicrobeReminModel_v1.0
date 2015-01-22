#Change working directory to MicrobeReminModel_v1.0/RCode using the setwd() command or using the Misc menu.  
#getwd() #shows the current working directory for the R gui.

filepath = paste("../Output/MRM0D/", sep="")

#Bacterial growth is assessed after 3 hours of intial spin-up.  Bacterial production is calculated over a 9 minute time interval which is the residence of exoenzyme in the particle.  The time step is 30 seconds. 
tsegs1<-362  #Begin time for bacterial production calculation
tsegs2<-379  #End time for bacterial production calculation

#--------------------------
#  Find Epsilon Optimum
#--------------------------
folderstate<-paste(filepath, "StateVariables/", sep="")
statenames<-list.files(folderstate)

folderbrates<-paste(filepath, "BacteriaRates/", sep="")
bratesnames<-list.files(folderbrates)

foldermrates<-paste(filepath, "MassTransferRates/", sep="")
mratesnames<-list.files(foldermrates)

analysis<-as.data.frame(matrix(NA, length(statenames), 6))
colnames(analysis)<-c("maxe", "qfrac.init", "ba.init", "ba.prod", "enz.prod", "h.prod")

for(a in 1:length(statenames)){

state<-read.table(paste(folderstate, statenames[a], sep=""))
colnames(state)<-c("time", "pom", "enz", "h", "ba", "hdom", "bf", "edom", "denz", "dedom")

brates<-read.table(paste(folderbrates, bratesnames[a], sep=""))
colnames(brates)<-c("time", "attach", "detach", "ba_grow", "bf_grow", "up_ba", "up_bf", "epsilon")

mrates<-read.table(paste(foldermrates, mratesnames[a], sep=""))
colnames(mrates)<-c("time", "pnum", "pom_dgd", "h_hdom", "ehl_k", "enz_dom", "denz_dom")

analysis$maxe[a]<-as.numeric(substr(statenames[a], 20, 25))
analysis$qfrac.init[a]<-as.numeric(substr(statenames[a], 32, 35))
analysis$ba.init[a]<-state$ba[1]
cperb<-5E-12  #carbon per bacteria
analysis$ba.init[a]<-analysis$ba.init[a]/cperb/mrates$pnum[2] #number of bacteria per particle

state<-state[tsegs1:tsegs2,]
brates<-brates[tsegs1:tsegs2,]
mrates<-mrates[tsegs1:tsegs2,]
analysis$ba.prod[a]<-sum(state$ba*brates$ba_grow)
analysis$enz.prod[a]<-sum(state$ba*brates$epsilon)
analysis$h.prod[a]<-sum(state$enz*mrates$pom_dgd)

}

analysis$qfrac.init<-analysis$qfrac.init*100

opt<-aggregate(analysis$ba.prod, list(analysis$qfrac.init), max)
colnames(opt)<-c("qfrac.init", "ba.prod")
opt<-merge(opt, analysis)
colnames(opt)<-gsub("maxe", "opte", colnames(opt))
opt<-opt[order(opt$qfrac.init),]

opt<-opt[c("ba.init", "qfrac.init", "opte")]

outfile<-paste("../Analysis/OptimalEpsilon.txt", sep="")
write.table(opt, outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


#---------------------------
#   Plot Optimal Epsilon
#---------------------------

OutGraph<-paste("../Graphs/MRM0D1_OptimalEpsilon.ps", sep="")
postscript(OutGraph, family="Helvetica", width=4, height=3.25)

#quartz(width=4, height=3.25)
par(oma=c(3,3,1,1))
par(mar=c(1.5,1,1,0.5)+0.1)
par(mfrow=c(1,1))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=1)
par(cex.lab=1)
par(cex.main=1)


plot(opt$ba.init, opt$opte, type="l", ylim=c(0,0.05), xaxt="n", yaxt="n", lwd=2)
axis(side=1, labels=TRUE)
axis(side=2, labels=TRUE, las=1)


mtext("Initial bacterial abundance per particle", side=1, cex=1, line=0.7 , adj=0.5, outer=TRUE)
mtext(expression(paste("Optimal  ", epsilon)), side=2, cex=1, line=1.6 , adj=0.5, outer=TRUE)


dev.off()

