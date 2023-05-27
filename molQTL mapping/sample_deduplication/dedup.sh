#!/bin/bash

set -e

### Load software
module load plink/1.90
module load bcftools/1.10.2

### parse params
tissue=$1


### Assign variant_id
if [ ! -f Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.vcf.gz ]
then
  bcftools annotate --threads 12 --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -a /group/zhougrp2/dguan/00_ref/gallus_gallus.vcf.gz -c ID -O z -o Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.vcf.gz /group/zhougrp2/dguan/23_snp_imputation/hjzhou98grp/Chicken_GTEx_RNAseq_SNPs_imputed_filtered.vcf.gz
  tabix -p vcf -f Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.vcf.gz
fi

### Convert vcf to plink format
if [ ! -f Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.bed ]
then
  plink --keep-allele-order --chr-set 33 --double-id --vcf Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.vcf.gz --make-bed --out Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated
fi
cat All_samples_mapping_rate.txt | awk '{if ($4 >= 1000000) print $1"\t"$1}' > samples_with_1M_reads.txt
if [ ! -f Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.filteredByReads.bed ]
then
  plink --bfile Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated --keep-allele-order --chr-set 33 --keep samples_with_1M_reads.txt --make-bed --out Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.filteredByReads
fi

### Split samples by tissue
mkdir -p ${tissue}
cat ${tissue}.list | awk '{print $1"\t"$1}' > ${tissue}/${tissue}.list
plink --bfile Chicken_GTEx_RNAseq_SNPs_imputed_filtered.IDformated.filteredByReads --keep-allele-order --chr-set 33 --keep ${tissue}/${tissue}.list --make-bed --out ${tissue}/${tissue}

### Calculate IBD values between samples
plink --bfile ${tissue}/${tissue} --keep-allele-order --chr-set 33 --genome --out ${tissue}/${tissue}.IBD

### Step1: filtering individuals that show high relatedness with multiple individuals by IBD value
n_dep=`cat ${tissue}/${tissue}.IBD.genome | awk '{if (NR >=2 && $12 >= 0.9) print}' | wc -l`
function dedup_samples {
  cat ${tissue}/${tissue}.IBD.genome | awk '{if (NR >=2 && $12 >= 0.9) print $1"\t"$3}' > ${tissue}/${tissue}.highIBD.pairs.txt
  n=`cat ${tissue}/${tissue}.highIBD.pairs.txt | wc -l`
  if [ -f ${tissue}/${tissue}.dup.txt ];
  then
    rm ${tissue}/${tissue}.dup.txt
  fi
  for i in $(seq 1 $n);
  do
    depS1=`cat ${tissue}/${tissue}.highIBD.pairs.txt | awk -v n=$i '{if (NR == n) print $1}'`
    depS2=`cat ${tissue}/${tissue}.highIBD.pairs.txt | awk -v n=$i '{if (NR == n) print $2}'`
    reads1=`cat All_samples_mapping_rate.txt | awk -v s=${depS1} '{if ($1 == s) print $4}'`
    reads2=`cat All_samples_mapping_rate.txt | awk -v s=${depS2} '{if ($1 == s) print $4}'`
    if [ "$reads1" -lt "$reads2" ];
    then
      echo -e "${depS1}" >> ${tissue}/${tissue}.dup.txt
    else
      echo -e "${depS2}" >> ${tissue}/${tissue}.dup.txt
    fi
  done
  if [ -s ${tissue}/${tissue}.dup.txt ]
  then
    cat ${tissue}/${tissue}.dup.txt | sort | uniq | awk '{print $1"\t"$1}' > ${tissue}/${tissue}.dup_removed.list
    rm ${tissue}/${tissue}.dup.txt
  fi
  plink --bfile ${tissue}/${tissue} --keep-allele-order --remove ${tissue}/${tissue}.dup_removed.list --chr-set 33 --make-bed --out ${tissue}/${tissue}.dedup
  mv ${tissue}/${tissue}.dedup.bed ${tissue}/${tissue}.bed
  mv ${tissue}/${tissue}.dedup.fam ${tissue}/${tissue}.fam
  mv ${tissue}/${tissue}.dedup.bim ${tissue}/${tissue}.bim
  plink --bfile ${tissue}/${tissue} --keep-allele-order --chr-set 33 --genome --out ${tissue}/${tissue}.IBD
  n_dep=`cat ${tissue}/${tissue}.IBD.genome | awk '{if (NR >=2 && $12 >= 0.9) print}' | wc -l`
  return ${n_dep}
}

while [ ${n_dep} -gt 0 ]
do 
  dedup_samples
done

cat ${tissue}/${tissue}.fam | awk '{print $1"\t"$1}' > ${tissue}/${tissue}.final.list
n_sample=`cat ${tissue}/${tissue}.final.list | wc -l`
if [ ${n_sample} -ge 40 ]
then
  mkdir -p /group/zhougrp2/dguan/90_updates/02_genotypes/${tissue}/
  echo -e "${tissue}\t${n_sample}" >> /group/zhougrp2/dguan/90_updates/Tissues_with_avail_samples_for_eqtl_mapping.txt
  plink --bfile ${tissue}/${tissue} --keep-allele-order --allow-extra-chr --chr-set 33 --keep ${tissue}/${tissue}.final.list --make-bed --out /group/zhougrp2/dguan/90_updates/02_genotypes/${tissue}/ChickenGTEx.${tissue}
  echo -e "\n\n\n*************\t${tissue} has ${n_sample} available for eQTL mapping!\n\n\n"
fi

