#!/bin/bash
#  365000132000018
for cluster_num in 360009933000018 373000139000019 365000132000018
do
    for file in /n03data/ellien/Euclid_ICL/simulations/out3/${cluster_num}/*
    do
      njobs=$((${njobs}+1))
      echo "Launch Dawis on file ${file}"
      n=$(basename "$file")
      ncl=${n:0:-5}
      indir=/n03data/ellien/Euclid_ICL/simulations/out3/${cluster_num}
      outdir=/n03data/ellien/Euclid_ICL/wavelets/out3/${cluster_num}/${ncl}

      bash slurm_dawis.sh $indir $n $outdir
      
   done
done
