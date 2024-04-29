#!/bin/bash
#
# Example SLURM job script for ALICE

#SBATCH --job-name=STAR_index
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=60gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ra500@student.le.ac.uk

# loading star version 2.7.10b

module load star/2.7.10b-m3zkpic

#running star to index the hg38 analysis set using hg38.refGenes.gtf
STAR --runMode genomeGenerate \ # genomeGenerate is the parameter to run the indexing
     --genomeDir /scratch/alice/r/ra500/seq \ # the directory output of the indexing files
     --genomeFastaFiles /scratch/alice/r/ra500/hg38/hg38.analysisSet.fa \ # the path to the reference genome file
     --sjdbGTFfile /scratch/alice/r/ra500/hg38/hg38.refGene.gtf \ # the path to the annotation file
     --genomeChrBinNbits 10 \ # to optimise memory usage during the indexing
