#!/bin/bash

for file in /n03data/ellien/Euclid_ICL/simulations/out4/*

do
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      qsub qsub_dawis.sh -v ncl=${n:3:3},nf=${n}
done
