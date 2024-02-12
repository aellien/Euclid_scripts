#!/bin/bash

for cluster_num in 54000066009352
do
    for file in /n03data/ellien/Euclid_ICL/simulations/out1/${cluster_num}/*
    do
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      ncl=${num_cluster}_${n:0:-5}
      qsub qsub_dawis.sh -v ncl=NIR,nf=${n}
      done
done
