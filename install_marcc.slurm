#!/bin/sh
#SBATCH -J "install_squid"
#SBATCH -o install_squid.o%j.log
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p shared
#SBATCH -t 12:0:0

module reset
module load python/2.7-anaconda

date

/software/apps/anaconda/5.2/python/2.7/bin/python -u /scratch/groups/pclancy3/programs/squid/install_marcc.py > /scratch/groups/pclancy3/programs/squid/install_marcc.log 2>&1

date
