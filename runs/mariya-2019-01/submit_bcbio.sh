#!/bin/sh
#SBATCH -p medium
#SBATCH -J bcbio_O2
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 4-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000

export PATH=/n/app/bcbio/tools/bin:$PATH
bcbio_nextgen.py ../config/alignment.yaml -n 24 -t ipython -s slurm -q medium -r t=4-00:00 --timeout 2000
