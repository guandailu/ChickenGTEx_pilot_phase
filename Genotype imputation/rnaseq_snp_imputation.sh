#!/bin/bash -l

# Aouthor: Dailu Guan
# Date: August 19, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=rnaseq_snp_imputation
#SBATCH --partition=bmm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=480G
#SBATCH --time=30-00:00:00
#SBATCH --output=/group/zhougrp2/dguan/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp2/dguan/98_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dguan@ucdavis.edu

module load bcftools/1.10.2
module load samtools/1.9
chr=$1


###### Keep final RNAseq samples
echo "Running RNAseq sample filtering"
bcftools view --threads 64 -i "F_MISSING<0.5" --min-af 0.05:minor -m2 -M2 -S /group/zhougrp2/dguan/15_RNAseq_filteredVCF/Final_samples.list -O z -o /group/zhougrp2/dguan/15_RNAseq_filteredVCF/chr${chr}/Chicken_GTEx_RNAseq_SNPs_chr${chr}.test.vcf.gz /group/zhougrp2/dguan/15_RNAseq_filteredVCF/chr${chr}/Chicken_GTEx_filtrated.${chr}.SNPs.vcf.gz

###### SNP imputation
echo "Running imputation step"
time java -Xmx480G -jar /home/dguan/bin/beagle.18May20.d20.jar gt=/group/zhougrp2/dguan/15_RNAseq_filteredVCF/chr${chr}/Chicken_GTEx_RNAseq_SNPs_chr${chr}.test.vcf.gz ref=/group/hjzhou98grp/dguan/40_WGS/10_wgs_PhasedSNPs/Chicken_GTEx_chr${chr}.phased.vcf.gz impute=true out=/home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed seed=666666 nthreads=64
time tabix -p vcf /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.vcf.gz


##### Filter SNPs with R2 > 0.8 and maf > 0.05
echo "Filtering phased SNPs"
time bcftools view -i 'DR2>0.8' --min-af 0.05:minor -l 9 --threads 64 -O z -o /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.filtered.vcf.gz /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.vcf.gz
time tabix -f -p vcf /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.filtered.vcf.gz
time md5sum /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.filtered.vcf.gz > /home/dguan/23_snp_imputation/Chicken_GTEx_RNAseq_SNPs_chr${chr}.imputed.filtered.vcf.gz.md5