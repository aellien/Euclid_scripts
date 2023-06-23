#!/bin/bash

for file in /n03data/ellien/Euclid_ICL/simulations/*crop*

do
      echo "Launch Dawis on file ${file}"
      echo "qsub qsub_dawis.sh -v ncl=${file:3:6},nf=${file}"
done
