#!/bin/bash -l
# Aouthor: Dailu Guan
# Date: August 19, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=Phase_reference_panel
#SBATCH --partition=bmm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=30-00:00:00
#SBATCH --output=/group/zhougrp2/dguan/40_WGS/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp2/dguan/40_WGS/98_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dguan@ucdavis.edu

module load bcftools/1.10.2
module load samtools/1.9
chr=$1

time bcftools view -m 2 -M 2 --min-af 0.05:minor -i 'F_MISSING<0.1' --threads 24 -l 9 -O z -o /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.vcf.gz /home/dguan/09_wgs_filteredSNPs/Chicken_GTEx_chr${chr}.final.vcf.gz

time zcat /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.vcf.gz | sed "s/\s\.:/\t.\/.:/g" > /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf
bgzip --threads 24 -l 9 /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf

time tabix -p vcf /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf.gz

####### PHASE REFERENCE PANEL
time java -Xmx120G -jar /home/dguan/bin/beagle.18May20.d20.jar gt=/home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf.gz out=/home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.phased seed=666666 nthreads=24
time tabix -p vcf /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.phased.vcf.gz

mv /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf.gz /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.MAFiltered.vcf.gz
mv /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.formated.vcf.gz.tbi /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.MAFiltered.vcf.gz.tbi
rm /home/dguan/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.biallelicSNPs.vcf.gz