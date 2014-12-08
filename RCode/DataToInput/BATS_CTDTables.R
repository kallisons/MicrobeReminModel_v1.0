#Change working directory to MicrobeReminModel_v1.0 using the setwd() command or using the Misc menu.  
getwd() #shows the current working directory for the R gui.  

#---------------------
# Source R Functions
#---------------------
source("RCode/R_functions/sigmat.fun")
source("RCode/R_functions/thetap.fun")
source("RCode/R_functions/alr.fun")

#------------------------------
# Read in List of Data Files
#------------------------------

DataFolder<-paste("Data/BATS/CTD/", sep="")
DataFiles<-list.files(DataFolder)

#-----------------
# Create Tables
#-----------------

depths<-seq(150, 4000, 10) # Final data in 10 m intervals

TempTable<-as.data.frame(matrix(NA, length(depths), 1))
colnames(TempTable)<-c("Depth")
TempTable$Depth<-depths
SalTable<-as.data.frame(matrix(NA, length(depths), 1))
colnames(SalTable)<-c("Depth")
SalTable$Depth<-depths
RhoTable<-as.data.frame(matrix(NA, length(depths), 1))
colnames(RhoTable)<-c("Depth")
RhoTable$Depth<-depths

#----------------------------
# Read and Process CTD Data
#----------------------------

for(k in 1:length(DataFiles)){

FileName<-paste(DataFolder, DataFiles[k], sep="")

CTD<-read.table(FileName)

colnames(CTD)<-c("ID", "Date", "Lat", "Lon", "Pres", "Depth", "Temp", "Sal", "DO", "Beam", "Fluor")

CTD<-CTD[c("ID", "Date", "Pres","Depth","Temp","Sal")]

#Replace Bad Values with NA
CTD$Sal<-ifelse(CTD$Sal<0, NA, CTD$Sal)
CTD$Sal<-ifelse(CTD$Sal>38, NA, CTD$Sal)
CTD$Temp<-ifelse(CTD$Temp<0, NA, CTD$Temp)

#Time Process
CTD$Year<-floor(CTD$Date)
CTD$Yfrac<-CTD$Date-CTD$Year
CTD$Sfrac<-unclass(ISOdate(CTD$Year+1,1,1,0,0,0)) - unclass(ISOdate(CTD$Year,1,1,0,0,0))
CTD$Date2<-ISOdate(CTD$Year,1,1) + CTD$Yfrac * CTD$Sfrac
CTD$Date<-substr(CTD$Date2,1,10)

CTD$Cast<-substr(CTD$ID, 6, 7)

#Average CTD casts done on the same day
CTDagg<-aggregate(CTD[c("Pres", "Temp", "Sal")], list(CTD[,c("Year")], CTD[,c("Date")], CTD[,c("Depth")]), mean, na.rm=TRUE)
colnames(CTDagg)<-c("Year","Date","Depth","Pres", "Temp", "Sal")
CTDagg$Depth<-round(CTDagg$Depth, digits=0)

#Interpolate CTD data
CTDinterp<-as.data.frame(matrix(NA, max(CTDagg$Depth), 2))
colnames(CTDinterp)<-c("Depth","Blank")
CTDinterp$Depth<-as.numeric(row.names(CTDinterp))
CTDinterp<-subset(CTDinterp, CTDinterp$Depth>=min(CTDagg$Depth, na.rm=TRUE))
CTDinterp2<-merge(CTDinterp, CTDagg, all=TRUE)
CTDinterp2<-CTDinterp2[,-2]
p1<-approx(CTDinterp2$Depth, CTDinterp2$Pres, CTDinterp2$Depth)
CTDinterp2$Pres<-p1$y
t1<-approx(CTDinterp2$Depth, CTDinterp2$Temp, CTDinterp2$Depth)
CTDinterp2$Temp<-t1$y
if(length(table(CTDinterp2$Sal))>0){
s1<-approx(CTDinterp2$Depth, CTDinterp2$Sal, CTDinterp2$Depth)
CTDinterp2$Sal<-s1$y
}

#Calculate density
CTDinterp2$Rho<-sigmat(CTDinterp2$Temp, CTDinterp2$Sal, CTDinterp2$Pres)

#Add time to the interpolations
CTDinterp2<-CTDinterp2[order(CTDinterp2$Pres),]
CTDinterp2$Year<-rep(CTDinterp2$Year[1])
CTDinterp2$Date<-rep(CTDinterp2$Date[1])
CTDinterp2<-unique(CTDinterp2)

#Move data to final tables
depths<-seq(150, 4000, 10)
depths<-subset(depths, depths<max(CTDinterp2$Depth))
diff<-(nrow(TempTable)-length(depths))
for(d in 1:length(depths)){
depths[d]<-which(CTDinterp2$Depth==eval(depths[d]))
}
CTDfinal<-CTDinterp2[c(depths),]
Date<-CTDfinal$Date[1]
TempTable<-cbind(TempTable, c(CTDfinal$Temp, rep(NA, eval(diff))))
colnames(TempTable)[k+1]<-eval(Date)
SalTable<-cbind(SalTable, c(CTDfinal$Sal, rep(NA, eval(diff))))
colnames(SalTable)[k+1]<-eval(Date)
RhoTable<-cbind(RhoTable, c(CTDfinal$Rho, rep(NA, eval(diff))))
colnames(RhoTable)[k+1]<-eval(Date)

}

#-------------
# Bad Data
#-------------
#Bad Salinity Values on 1989-06, replace with NA
SalTable[,10]<-rep(NA)
RhoTable[,10]<-rep(NA)
#Bad Salinity Values on 1990-03, replace with NA
SalTable[,18]<-rep(NA)
RhoTable[,18]<-rep(NA)

#----------------------------
# Plot Profiles of CTD Data
#----------------------------
quartz()
plot(TempTable[,2], TempTable$Depth*-1, type="l", xlim=c(0,22))
last<-ncol(TempTable)-2
for(y in 3:last){
	lines(TempTable[,y], TempTable$Depth*-1)}

quartz()
plot(SalTable[,2], SalTable$Depth*-1, type="l")
last<-ncol(SalTable)-2
for(y in 3:last){
	lines(SalTable[,y], SalTable$Depth*-1)}

quartz()
plot(RhoTable[,2], RhoTable$Depth*-1, type="l")
last<-ncol(RhoTable)-2
for(y in 3:last){
	lines(RhoTable[,y], RhoTable$Depth*-1)}

#------------------------------------
# Calculate Long-term Means of Data
#------------------------------------
TempTable2<-TempTable
SalTable2<-SalTable
RhoTable2<-RhoTable

TempTable$mean<-rowMeans(TempTable[,2:ncol(TempTable)], na.rm=TRUE)

TempTable[,2:length(TempTable)]<-round(TempTable[,2:length(TempTable)], digits=1)

SalTable$mean<-rowMeans(SalTable[,2:ncol(SalTable)], na.rm=TRUE)

SalTable[,2:9]<-round(SalTable[,2:9], digits=1)
SalTable[,11:17]<-round(SalTable[,11:17], digits=1)
SalTable[,19:length(SalTable)]<-round(SalTable[,19:length(SalTable)], digits=1)

RhoTable$mean<-rowMeans(RhoTable[,2:ncol(RhoTable)], na.rm=TRUE)
RhoTable[,2:9]<-round(RhoTable[,2:9], digits=2)
RhoTable[,11:17]<-round(RhoTable[,11:17], digits=2)
RhoTable[,19:length(RhoTable)]<-round(RhoTable[,19:length(RhoTable)], digits=2)

#--------------------------------
# Write Processed Data to Files
#--------------------------------
outfile<-paste("Data/BATS/BATS_CTDTemp_10m.txt", sep="")
write.table(TempTable, outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=" ")

outfile<-paste("Data/BATS/BATS_CTDSal_10m.txt", sep="")
write.table(SalTable, outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=" ")

outfile<-paste("Data/BATS/BATS_CTDRho_10m.txt", sep="")
write.table(RhoTable, outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=" ")

