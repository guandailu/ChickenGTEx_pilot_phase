#!/bin/bash


prefix=$1
chr=$2
sampleinfo=$3
#### Activate software
source ~/.bashrc
conda activate eigensoft
module load plink

#### Convert vcf to plink
#plink --vcf ${prefix}.vcf.gz --chr-set ${chr} --double-id --recode --out ${prefix}


#### config a file for converting plink format to EIGENSTRAT format
#genotypename:    Chicken_GTEx_RNAseq_SNPs_imputed_filtered.ped
#snpname:         Chicken_GTEx_RNAseq_SNPs_imputed_filtered.map # or example.map, either works 
#indivname:       Chicken_GTEx_RNAseq_SNPs_imputed_filtered.ped # or example.ped, either works
#outputformat:    EIGENSTRAT
#genotypeoutname: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.eigenstratgeno
#snpoutname:      Chicken_GTEx_RNAseq_SNPs_imputed_filtered.snp
#indivoutname:    Chicken_GTEx_RNAseq_SNPs_imputed_filtered.ind
#familynames:     NO
#echo -e "genotypename: ${prefix}.ped" >> ${prefix}.ped2eigenstrat.par
#echo -e "snpname: ${prefix}.map" >> ${prefix}.ped2eigenstrat.par
#echo -e "indivname: ${prefix}.ped" >> ${prefix}.ped2eigenstrat.par
#echo -e "outputformat: EIGENSTRAT" >> ${prefix}.ped2eigenstrat.par
#echo -e "genotypeoutname: ${prefix}.eigenstratgeno" >> ${prefix}.ped2eigenstrat.par
#echo -e "snpoutname: ${prefix}.snp" >> ${prefix}.ped2eigenstrat.par
#echo -e "indivoutname: ${prefix}.ind" >> ${prefix}.ped2eigenstrat.par
#echo -e "familynames: NO" >> ${prefix}.ped2eigenstrat.par


#### run program for converting plink format to EIGENSTRAT format
#convertf -p ${prefix}.ped2eigenstrat.par

### config a file for running smartpca
#genotypename: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.geno
#snpname: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.snp
#indivname: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.evec
#evecoutname: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.pca.evec
#evaloutname: Chicken_GTEx_RNAseq_SNPs_imputed_filtered.eval
#altnormstyle: NO
#numoutevec: 2
#numoutlieriter: 5
#numoutlierevec: 2
#outliersigmathresh: 6.0
#qtmode: 0
#echo -e "genotypename: ${prefix}.eigenstratgeno" >> ${prefix}.pca.par
#echo -e "snpname: ${prefix}.snp" >> ${prefix}.pca.par
#echo -e "indivname: ${prefix}.ind" >> ${prefix}.pca.par
#echo -e "evecoutname: ${prefix}.evec" >> ${prefix}.pca.par
#echo -e "evaloutname: ${prefix}.eval" >> ${prefix}.pca.par
#echo -e "altnormstyle: NO" >> ${prefix}.pca.par
#echo -e "numoutevec: 50" >> ${prefix}.pca.par
#echo -e "numoutlieriter: 0" >> ${prefix}.pca.par
#echo -e "numoutlierevec: 10" >> ${prefix}.pca.par
#echo -e "outliersigmathresh: 6.0" >> ${prefix}.pca.par
#echo -e "qtmode: 0" >> ${prefix}.pca.par

### running smartpca program
#smartpca -p ${prefix}.pca.par > ${prefix}.pca.log

### Running plot script
Rscript smartpca_plot.R ${prefix} ${sampleinfo}
