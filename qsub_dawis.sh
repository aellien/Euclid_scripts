#!/bin/bash
#PBS -o /home/ellien/Euclid_ICL/logs/${ncl}.out
#PBS -j oe
#PBS -N icl_euclid
#PBS -l nodes=1:ppn=4,walltime=47:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
echo ${ncl}
echo ${nf}
python /home/ellien/Euclid_ICL/Euclid_scripts/dawis_${ncl}.py ${nf}

exit 0
