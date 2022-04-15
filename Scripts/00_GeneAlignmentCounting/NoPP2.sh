#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J NoPP2
#SBATCH -o OUT/NoPP2.out
#SBATCH -e ERR/NoPP2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=NoPP2 \
                   --transcriptome=../PIGscRNAseq1/ssc97 \
                   --fastqs=Nopp2 \
                   --sample=NoPP2-1,NoPP2-2,NoPP2-3,NoPP2-4 \
                   --localcores=38

