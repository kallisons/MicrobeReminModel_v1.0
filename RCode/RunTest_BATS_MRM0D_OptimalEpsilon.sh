#!/bin/sh

echo “——————————————————————————————”

echo “MRM0D Optimal Epsilon”

diff ../Analysis/OptimalEpsilon.txt ../TestFiles/MRM0D/OptimalEpsilon/OptimalEpsilon.txt
if [ $? -eq 0 ]
  then
    echo “** PASSES TEST: Output Files Are The Same **”
  else
    echo “** FAILS TEST: Output Files Are Different **”
fi

echo “——————————————————————————————”