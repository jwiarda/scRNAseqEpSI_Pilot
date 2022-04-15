#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J JEJ2
#SBATCH -o OUT/JEJ2.out
#SBATCH -e ERR/JEJ2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=JEJ2 \
                   --transcriptome=../PIGscRNAseq1/ssc97 \
                   --fastqs=Jej2 \
                   --sample=Jej2-1,Jej2-2,Jej2-3,Jej2-4 \
                   --localcores=38

