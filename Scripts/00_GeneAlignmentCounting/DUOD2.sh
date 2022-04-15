#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J DUOD2
#SBATCH -o OUT/DUOD2.out
#SBATCH -e ERR/DUOD2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=DUOD2 \
                   --transcriptome=../PIGscRNAseq1/ssc97 \
                   --fastqs=Duod2 \
                   --sample=Duod2-1,Duod2-2,Duod2-3,Duod2-4 \
                   --localcores=38

