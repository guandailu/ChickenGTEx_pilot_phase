#!/bin/bash



gwas_file=$1
tissue=$2
qtl_type=$3


metaimputation_gwas_files=("RIRline.E3.assoc.mlma.fdr.txt" "RIRline.BSPO2.assoc.mlma.fdr.txt" "RIRline.EM.assoc.mlma.fdr.txt" "RIRline.EW1.assoc.mlma.fdr.txt" "WLline.BW4.assoc.mlma.fdr.txt" "WLline.CVEW.assoc.mlma.fdr.txt" "RIRline.SP1.assoc.mlma.fdr.txt" "WLline.YW1.assoc.mlma.fdr.txt" "WLline.AFC.assoc.mlma.fdr.txt")

stepwise_gwas_files=("WLline.CO1.assoc.mlma.fdr.txt" "WLline.EW1.assoc.mlma.fdr.txt" "WLline.AH1.assoc.mlma.fdr.txt" "RIRline.BW4.assoc.mlma.fdr.txt" "RIRline.AH2.assoc.mlma.fdr.txt" "WLline.AH2.assoc.mlma.fdr.txt" "RIRline.AFC.assoc.mlma.fdr.txt" "WLline.E3.assoc.mlma.fdr.txt" "RIRline.AH1.assoc.mlma.fdr.txt" "WLline.CO2.assoc.mlma.fdr.txt" "WLline.SP1.assoc.mlma.fdr.txt" "WLline.YW2.assoc.mlma.fdr.txt" "WLline.EW2.assoc.mlma.fdr.txt" "RIRline.BW2.assoc.mlma.fdr.txt" "RIRline.GEGG2.assoc.mlma.fdr.txt" "WLline.BW2.assoc.mlma.fdr.txt" "WLline.BSPO1.assoc.mlma.fdr.txt" "RIRline.GEGG3.assoc.mlma.fdr.txt" "WLline.BSE1.assoc.mlma.fdr.txt" "WLline.AFE.assoc.mlma.fdr.txt" "RIRline.YW2.assoc.mlma.fdr.txt" "WLline.DF1.assoc.mlma.fdr.txt" "RIRline.EW2.assoc.mlma.fdr.txt" "WLline.GEGG2.assoc.mlma.fdr.txt" "RIRline.BSE1.assoc.mlma.fdr.txt" "RIRline.SP2.assoc.mlma.fdr.txt" "WLline.GEGG3.assoc.mlma.fdr.txt")

population="${gwas_file%%.*}"
[[ ${metaimputation_gwas_files[*]} =~ ${gwas_file} ]] && genotype_file="../GWAS_results/Meta-imputation/Imputed_Results/${population}.MetaImputed.IDformatted.filtered" || genotype_file="../GWAS_results/Stepwise_imputation/Imputed_variants/${population}.imputed.filtered.IDformatted"

if [[ ${qtl_type} == "eQTL" ]]
then
    output_path="fastenloc_eqtl_res"
    fastenloc_ann="fastenloc_ann/coloc_eqtl/${tissue}.fastenloc.eqtl.annotation.vcf.gz"
elif [[ ${qtl_type} == "sQTL" ]]
then
    output_path="fastenloc_sqtl_res"
    fastenloc_ann="fastenloc_ann/coloc_sqtl/${tissue}.fastenloc.eqtl.annotation.vcf.gz"
elif [[ ${qtl_type} == "lncQTL" ]]
then
    output_path="fastenloc_lncqtl_res"
    fastenloc_ann="fastenloc_ann/coloc_lncqtl/${tissue}.fastenloc.eqtl.annotation.vcf.gz"
elif [[ ${qtl_type} == "exQTL" ]]
then
    output_path="fastenloc_exqtl_res"
    fastenloc_ann="fastenloc_ann/coloc_exonqtl/${tissue}.fastenloc.eqtl.annotation.vcf.gz"
elif [[ ${qtl_type} == "3aQTL" ]]
then
    output_path="fastenloc_3aqtl_res"
    fastenloc_ann="fastenloc_ann/coloc_3aqtl/${tissue}.fastenloc.eqtl.annotation.vcf.gz"
fi
mkdir -p ${output_path}


for cmd in CalcZscore CalcPIP fastenloc;
do
   Rscript run_fastenloc.R ${cmd} -f ${genotype_file} -g ../Selected_GWAS/${gwas_file} --snp SNP --is_zscore FALSE --beta b --se se --plink /share/apps/plink-1.90/plink --torus /home/dguan/bin/torus/src/torus.static --fastenloc /share/apps/fastenloc-2.0/bin/fastenloc --fastenloc_ann ${fastenloc_ann} --tissue ${tissue} --thread 64 --out_prefix ${output_path}/${gwas_file%.assoc.mlma.fdr.txt*}
done
