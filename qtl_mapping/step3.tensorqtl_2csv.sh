#!/bin/bash

# Author: Dailu Guan
# Date: October 11, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=tensor
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/group/zhougrp2/dguan/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp2/dguan/98_logs/%x-%j.err
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dguan@ucdavis.edu


source ~/.bashrc
set -e

# Parse parameters
while getopts t:n:m: option
do
    case "${option}"
        in
        t)tissue=${OPTARG};;
	n)node=${OPTARG};; # node used for running the command
    esac
done

# Create Log directory
mkdir -p Logs
mkdir -p Temp

# Run by chromosome
for chr in {1..28} {30..33}; do
  echo -e '#!/bin/bash\n' "python3 step3.tensorqtl_2csv.py ${tissue}/ChickenGTEx.${tissue}.cis_qtl_pairs.${chr}.parquet ${tissue}/ChickenGTEx.${tissue}.cis_qtl_pairs.${chr}.parquet.csv.gz" > Temp/step3.tensorqtl_2csv.${tissue}.${chr}.sh
sbatch -p ${node} -c 1 --mem=2G -t 1-0 --job-name="parquet2csv" -o Logs/step3.tensorqtl_2csv.${tissue}.chr${chr}.%j.out Temp/step3.tensorqtl_2csv.${tissue}.${chr}.sh
done

echo "Done!"
