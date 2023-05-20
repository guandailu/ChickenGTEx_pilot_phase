#!/bin/bash -l

# Author: Dailu Guan
# Date: October 11, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=gpu_tensor
#SBATCH --account=gpul 
#SBATCH --partition=gpul
#SBATCH --gres=gpu:titan:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --output=/group/zhougrp2/dguan/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp2/dguan/98_logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dguan@ucdavis.edu


#############################################
# ENV ACTIVATE
source ~/.bashrc
conda activate R4

# Parse parameters
while getopts g:p:c:t:o: option
do 
    case "${option}"
        in
        g)genotype=${OPTARG};;
        p)phenotype=${OPTARG};;
        c)covariates=${OPTARG};;
	t)tissue=${OPTARG};;
        o)outpath=${OPTARG}
    esac
done

## Creat output directory
mkdir -p ${outpath}

prefix=${outpath}/ChickenGTEx.${tissue}
# 1) cis-QTL mapping: compute cis nominal associations for all variant-gene pairs
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis_nominal \
            --covariates ${covariates} \
            --maf_threshold 0.05

# 2) cis-QTL mapping: compute cis permutations and define eGenes
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis \
            --covariates ${covariates} \
            --maf_threshold 0.05

# 3) perform fdr test
Rscript step1.tensorqtl_fdr.R ${tissue}/ChickenGTEx.${tissue}.cis_qtl.txt.gz ${tissue}/ChickenGTEx.${tissue}.cis_qtl.txt 0.05
cat ${tissue}/ChickenGTEx.${tissue}.cis_qtl.txt | awk '{if (NR==1 || $19 =="TRUE") print}'  > ${tissue}/ChickenGTEx.${tissue}.cis_qtl.fdr.txt

echo "Done!"
