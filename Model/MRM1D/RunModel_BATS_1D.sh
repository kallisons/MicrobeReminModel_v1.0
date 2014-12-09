#!/bin/sh

rm filename1.txt
rm filename2.txt

filepath1=../../Input/BATS/

	echo "${filepath1}Force_30sTimeStep.in">filename1.txt
	echo "${filepath1}Profile_Mean10m.in">filename2.txt

#Right now qfrac is modified in the fortran file
qfrac="0.25"

#w is sinking speed in m/h and r_top is the particle radius in m at 150 m
for w in 2.00
do
for r_top in 0.0010
do
	rm MassTransferRates.out
	rm StateVariables.out
	rm BacteriaRates.out
	rm sizesink.input

	echo "&control" > sizesink.input
	echo "w =" ${w} >> sizesink.input
	echo "r_top =" ${r_top} >> sizesink.input
	echo "/" >> sizesink.input
	echo " " >> sizesink.input

./MRM1D1 < sizesink.input
cp StateVariables.out ../../Output/MRM1D/StateVariables/StateVariables_Sink${w}_Size${r_top}_Qfrac${qfrac}.out
cp BacteriaRates.out ../../Output/MRM1D/BacteriaRates/BacteriaRates_Sink${w}_Size${r_top}_Qfrac${qfrac}.out
cp MassTransferRates.out ../../Output/MRM1D/MassTransferRates/MassTransferRates_Sink${w}_Size${r_top}_Qfrac${qfrac}.out
done
done