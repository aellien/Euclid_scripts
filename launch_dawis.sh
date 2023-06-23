#!/bin/bash

for file in /n03data/ellien/Euclid_ICL/simulations/*crop*

do
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      echo "qsub qsub_dawis.sh -v ncl=${n:3:6},nf=${file}"
done
