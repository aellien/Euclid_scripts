#!/bin/bash

for cluster_num in 365000132000018 360009933000018 373000139000019
do
    for file in /n03data/ellien/Euclid_ICL/simulations/out2/${cluster_num}/*
    do
      njobs=$((${njobs}+1))
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      ncl=${cluster_num}_${n:0:-5}
      indir=/n03data/ellien/Euclid_ICL/simulations/out2/${cluster_num}/
      outdir=/n03data/ellien/Euclid_ICL/wavelets/out2/${cluster_num}
      #qsub qsub_dawis.sh -v ncl=${ncl},nf=${n},indir=${indir},outdir=${outdir}
      bash slurm_dawis.sh $indir $file $outdir
      
   done
done
