#!/bin/sh

#max_epsilon is the maximum exoenzyme production rate

#qfrac is the fraction of the surface area of the spherical particle covered by a single layer of bacteria, values greater than 1 mean that there is a second layer of bacteria

for max_epsilon in 0.0000 0.0100 0.0200 0.0300 0.0400 0.0500
do
for qfrac in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.10
do
	rm MassTransferRates.out
	rm StateVariables.out
	rm BacteriaRates.out
	rm grow2.input
	echo "&control" > grow2.input
	echo "qfrac =" ${qfrac} >> grow2.input
	echo "max_epsilon =" ${max_epsilon} >> grow2.input
	echo "/" >> grow2.input
	echo " " >> grow2.input

./MRM0D1 < grow2.input
cp StateVariables.out ../../Output/MRM0D/StateVariables/StateVariables_maxe${max_epsilon}_qfrac${qfrac}.out
cp BacteriaRates.out ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe${max_epsilon}_qfrac${qfrac}.out
cp MassTransferRates.out ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe${max_epsilon}_qfrac${qfrac}.out
done
done