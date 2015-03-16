#Change working directory to MicrobeReminModel_v1.0/Model/MRM0D using the setwd() command or using the Misc menu.  
#getwd() #shows the current working directory for the R gui.  

#The code in this file compares the output produced by the 0-dimensional version of the Microbial Remineralization Model with test files to make sure that the output is the same.  If the files have the same data (to 3 significant digits), then the tests are passed.  A passed test means that the model is working!  If the tests fail, the model is not working and should be compiled again.

decimalplaces <- function(x) {      
       hold<-strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)
       hold2<-ifelse(hold=="character(0)", list(c("","")), hold)
       df<-data.frame(matrix(unlist(hold2), nrow=length(hold2), byrow=T))
       nchar(as.character(df[,2])) 
   };


TableCompare <- function(filename.new, filename.old) {

	#filename.new<-paste("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out")
	#filename.old<-paste("../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out")

	new <- read.table(filename.new)
	old <- read.table(filename.old)

#--------------------------------------------------------
# Round to 3 significant digits
#--------------------------------------------------------	
	new <- signif(new, digits=3)
	old <- signif(old, digits=3)

	#-----------------------------------------------------
	# Check if there are the same number of rows in files
	#-----------------------------------------------------
	rowcomp <- abs(nrow(new) - nrow(old))

	if (rowcomp > 0) {
		cat("***ERROR***\n")
		cat("*Different Numbers of Rows*\n")
		cat(paste(filename.new, "\n", sep = ""))
		cat(paste(filename.old, "\n", sep = ""))
		cat("--------------\n")
		stop("FAILS TEST, Exiting RScript")
	}

	#--------------------------------------------------------
	# Check if there are the same number of columns in files
	#--------------------------------------------------------
	colcomp <- abs(ncol(new) - ncol(old))

	if (colcomp > 0) {
		cat("***ERROR***\n")
		cat("*Different Numbers of Columns*\n")
		cat(paste(filename.new, "\n", sep = ""))
		cat(paste(filename.old, "\n", sep = ""))
		cat("--------------\n")
		stop("FAILS TEST, Exiting RScript")
	}

	#-----------------------------------------------------
	# Check values in tables
	#-----------------------------------------------------
	numcols <- ncol(new)

	for (i in 1:ncol(new)) {

		diffcols <- abs(new[, i] - old[, i])

		finddiffs <- which(diffcols > 0)

		#---------------------------------
		# Find rounding errors 
		#---------------------------------

		if (length(finddiffs > 0)) {
			vals <- new[finddiffs, i]

			thds <- c()

			for (k in 1:length(vals)) {
					thds <- c(thds, 1 * 10^(-decimalplaces(vals[k])))
			}

			diffvals <- diffcols[finddiffs]
			bigdiffs <- subset(finddiffs, diffvals > thds)

		#-----------------------------
		# Print non-rounding errors 
		#-----------------------------
			if (length(bigdiffs) > 0) {
				for (j in 1:length(bigdiffs)) {
					cat("***ERROR***\n")
					cat("*Different Values*\n")
					cat(paste("col=", i, "\t row=", bigdiffs[j], "\n", sep = ""))
					cat(paste("File 1:\n", sep = ""))
					cat(paste(filename.new, "\n", sep = ""))
					cat(paste("value=", new[bigdiffs[j], i], "\n", sep = ""))
					cat(paste("File 2:\n", sep = ""))
					cat(paste(filename.old, "\n", sep = ""))
					cat(paste("value=", old[bigdiffs[j], i], "\n", sep = ""))
					cat("--------------\n")
					stop("FAILS TEST, Exiting RScript")
				}
			}
		}
	}
	cat("***PASSES ALL COMPARISON TESTS***\n")
	cat(paste("File 1:\n", sep = ""))
	cat(paste(filename.new, "\n", sep = ""))
	cat(paste("File 2:\n", sep = ""))
	cat(paste(filename.old, "\n", sep = ""))
	cat("--------------\n")

}

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0100_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0100_qfrac0.50.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.50.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.75.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.75.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac1.00.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac1.00.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0300_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0300_qfrac0.50.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0400_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0400_qfrac0.50.out")

TableCompare("../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0500_qfrac0.50.out", "../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0500_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0000_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0000_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0100_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0100_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.75.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.75.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac1.00.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac1.00.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0300_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0300_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0400_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0400_qfrac0.50.out")

TableCompare("../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0500_qfrac0.50.out", "../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0500_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0000_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0000_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0100_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0100_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.75.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.75.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac1.00.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac1.00.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0300_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0300_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0400_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0400_qfrac0.50.out")

TableCompare("../../Output/MRM0D/StateVariables/StateVariables_maxe0.0500_qfrac0.50.out", "../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0500_qfrac0.50.out")
