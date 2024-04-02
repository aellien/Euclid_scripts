#!/bin/bash
#PBS -o /n03data/ellien/Euclid_ICL/logs/${ncl}.out
#PBS -j oe
#PBS -N icl_euclid
#PBS -l nodes=1:ppn=1,walltime=47:59:00
#PSB -S /bin/bash

#module load intelpython/3-2020.4
conda init bash
source /home/ellien/.bashrc
conda activate dawis

echo ${ncl}
echo ${nf}
python /home/ellien/Euclid_ICL/Euclid_scripts/dawis_NIR.py ${indir} ${nf} ${outdir}

exit 0
