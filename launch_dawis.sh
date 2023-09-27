#!/bin/bash

for file in /n03data/ellien/Euclid_ICL/simulations/out6/*bgsub*

do
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      #ncl=${n:3:3}
      qsub qsub_dawis.sh -v ncl=NIR,nf=${n}
done
