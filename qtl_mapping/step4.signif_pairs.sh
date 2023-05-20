#!/bin/bash

# Author: Dailu Guan
# Date: October 11, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=tensor
#SBATCH --partition=bmh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --output=/group/zhougrp2/dguan/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp2/dguan/98_logs/%x-%j.err
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dguan@ucdavis.edu
set -e

# Parse parameters
while getopts t: option
do
    case "${option}"
        in
        t)tissue=${OPTARG}
    esac
done

# Create Log directory
Rscript step4.signif_pairs.R ${tissue}

