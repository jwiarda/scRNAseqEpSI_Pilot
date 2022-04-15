#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J IPP2
#SBATCH -o OUT/IPP2.out
#SBATCH -e ERR/IPP2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=IPP2 \
                   --transcriptome=../PIGscRNAseq1/ssc97 \
                   --fastqs=Ipp2 \
                   --sample=IPP2-1,IPP2-2,IPP2-3,IPP2-4 \
                   --localcores=38

