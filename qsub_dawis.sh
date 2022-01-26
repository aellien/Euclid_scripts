#!/bin/bash
#PBS -o /home/ellien/Euclid_ICL/logs/icl_euclid.out
#PBS -j oe
#PBS -N icl_euclid
#PBS -l nodes=1:ppn=6,walltime=120:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
python /home/ellien/Euclid_ICL/Euclid_scripts/dawis_euclid.py


exit 0
