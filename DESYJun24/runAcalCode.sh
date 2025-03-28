#!/bin/bash
BOARD=$1

for C in {0..15}; do
   for R in {0..15}; do
       echo "python3 mainRDF.py --board $BOARD --col $C --row $R > logsAcalCode/log_B${BOARD}C${C}R${R}.log"
       python3 mainRDF.py --board $BOARD --col $C --row $R > logsAcalCode/log_B${BOARD}C${C}R${R}.log;
   done;
done
