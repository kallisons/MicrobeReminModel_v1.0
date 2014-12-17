#!/bin/sh

rm filename1.txt
rm filename2.txt

filepath1=../../Input/BATS/

	echo "${filepath1}Force_30sTimeStep.in">filename1.txt
	echo "${filepath1}Profile_Mean10m.in">filename2.txt

#w is sinking speed in m h^-1 at 150 m
#r_top is the particle radius in m at 150 m
#qfrac is fraction of the particle surface covered by attached bacteria at 150 m

for w in 2.00
do
for r_top in 0.0010
do
for qfrac in 0.50
do
	rm sizesink.input	

	echo "&control" > sizesink.input
	echo "w =" ${w} >> sizesink.input
	echo "r_top =" ${r_top} >> sizesink.input
	echo "qfrac =" ${qfrac} >> sizesink.input
	echo "/" >> sizesink.input
	echo " " >> sizesink.input

./MRM1D1 < sizesink.input
cp StateVariables.out ../../Output/MRM1D/StateVariables/StateVariables_Sink${w}_Size${r_top}_Qfrac${qfrac}.out
cp BacteriaRates.out ../../Output/MRM1D/BacteriaRates/BacteriaRates_Sink${w}_Size${r_top}_Qfrac${qfrac}.out
cp MassTransferRates.out ../../Output/MRM1D/MassTransferRates/MassTransferRates_Sink${w}_Size${r_top}_Qfrac${qfrac}.out

rm MassTransferRates.out
rm StateVariables.out
rm BacteriaRates.out

done
done
done

rm filename1.txt
rm filename2.txt
rm sizesink.input