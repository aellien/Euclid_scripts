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

ncpus=1

export OMP_NUM_THREADS=${ncpus}
export OPENBLAS_NUM_THREADS=${ncpus} 
export MKL_NUM_THREADS=${ncpus}
export VECLIB_MAXIMUM_THREADS=${ncpus}
export NUMEXPR_NUM_THREADS=${ncpus}

echo ${ncl}
echo ${nf}

python /home/ellien/Euclid_ICL/Euclid_scripts/dawis_NIR.py ${indir} ${nf} ${outdir}

exit 0