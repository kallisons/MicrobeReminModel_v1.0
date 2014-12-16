#!/bin/sh

echo “——————————————————————————————”

echo “Bacterial Rates”

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0000_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0100_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0100_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.75.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac0.75.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac1.00.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0200_qfrac1.00.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0300_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0300_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0400_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0400_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/BacteriaRates/BacteriaRates_maxe0.0500_qfrac0.50.out ../../TestFiles/MRM0D/BacteriaRates/BacteriaRates_maxe0.0500_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

echo “——————————————————————————————”

echo “Mass Transfer Rates”

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0000_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0000_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0100_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0100_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.75.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac0.75.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac1.00.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0200_qfrac1.00.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0300_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0300_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0400_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0400_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/MassTransferRates/MassTransferRates_maxe0.0500_qfrac0.50.out ../../TestFiles/MRM0D/MassTransferRates/MassTransferRates_maxe0.0500_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi


echo “——————————————————————————————”
echo “State Variables”

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0000_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0000_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0100_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0100_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.75.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac0.75.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac1.00.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0200_qfrac1.00.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0300_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0300_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0400_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0400_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

diff ../../Output/MRM0D/StateVariables/StateVariables_maxe0.0500_qfrac0.50.out ../../TestFiles/MRM0D/StateVariables/StateVariables_maxe0.0500_qfrac0.50.out
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi


echo “——————————————————————————————”