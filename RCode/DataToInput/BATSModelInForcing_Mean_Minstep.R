TrapC<-read.table("/Users/kasmith/Documents/EnzymeBacteriaModel/Data/BATS/BATS_trapC.txt", header=TRUE)
TrapC<-subset(TrapC, TrapC$yymmdd1>=19900000)

TrapC150<-subset(TrapC, TrapC$Dep==150)
TrapC150$month<-substr(TrapC150$yymmdd1, 5,6)
TrapC150$Cavg<-ifelse(TrapC150$Cavg==-9.99, NA, TrapC150$Cavg)

TrapC200<-subset(TrapC, TrapC$Dep==200)
TrapC200$month<-substr(TrapC200$yymmdd1, 5,6)
TrapC200$Cavg<-ifelse(TrapC200$Cavg==-9.99, NA, TrapC200$Cavg)

C150flux<-round(mean(TrapC150$Cavg, na.rm=TRUE), digits=1)/24/60

days<-365*1

makeinput<-as.data.frame(matrix(NA,days*24*60, 3))
colnames(makeinput)<-c("jday", "hhmm", "POC")

n<-1

for(x in 1:days){
p<-n+24*60-1

makeinput$jday[n:p]<-rep(x)
makeinput$hhmm[n:p]<-c(seq(0,59,1), seq(100,159,1), seq(200,259,1), seq(300,359,1), seq(400,459,1), seq(500,559,1), seq(600,659,1), seq(700,759,1), seq(800,859,1), seq(900,959,1), seq(1000,1059,1), seq(1100,1159,1), seq(1200,1259,1), seq(1300,1359,1), seq(1400,1459,1), seq(1500,1559,1), seq(1600,1659,1), seq(1700,1759,1), seq(1800,1859,1), seq(1900,1959,1), seq(2000,2059), seq(2100,2159,1), seq(2200,2259,1), seq(2300,2359,1))
n<-n+24*60

makeinput$POC<-round(rep(C150flux), digits=3)

outfile<-paste("/Users/kasmith/Documents/EnzymeBacteriaModel/BacteriaModel/Experiments/TrapData/BATSin/POM_In_MinStep.in", sep="")
write.table(makeinput, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
}




