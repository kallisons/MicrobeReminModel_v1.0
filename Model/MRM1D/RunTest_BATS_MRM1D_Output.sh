#!/bin/sh

echo “——————————————————————————————”

echo “Bacterial Rates”

diff ../../Output/MRM1D/BacteriaRates/BacteriaRates_Sink2.00_Size0.0010_Qfrac0.50.out ../../TestFiles/MRM1D/BacteriaRates/BacteriaRates_Sink2.00_Size0.0010_Qfrac0.25.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

echo “——————————————————————————————”

echo “Mass Transfer Rates”

diff ../../Output/MRM1D/MassTransferRates/MassTransferRates_Sink2.00_Size0.0010_Qfrac0.50.out ../../TestFiles/MRM1D/MassTransferRates/MassTransferRates_Sink2.00_Size0.0010_Qfrac0.25.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

echo “——————————————————————————————”
echo “State Variables”

diff ../../Output/MRM1D/StateVariables/StateVariables_Sink2.00_Size0.0010_Qfrac0.50.out ../../TestFiles/MRM1D/StateVariables/StateVariables_Sink2.00_Size0.0010_Qfrac0.25.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

echo “——————————————————————————————”