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
#SBATCH --time=4:00:00
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

prefix=${outpath}/ChickenGTEx.${tissue}
cis_output=${outpath}/ChickenGTEx.${tissue}.cis_qtl.fdr.txt
# 4) cis-QTL mapping: conditionally independent QTLs
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis_independent \
            --covariates ${covariates} \
            --cis_output ${cis_output} \
            --maf_threshold 0.05
# 5) trans-QTL mapping
python3 -m tensorqtl \
           ${genotype} ${phenotype} ${prefix} \
           --covariates ${covariates} \
           --mode trans \
           --maf_threshold 0.05 --output_text 

echo "Done!"
