#!/bin/bash
source ~/.bashrc

dir_output="./"
dir_nominal="/group/zhougrp2/dguan/90_updates/02_eQTLs/04_results/"
subset_size=1000000

tis_names=("Adipose" "Blood" "Brain" "Bursa" "Cecum" "Cerebellum" "Duodenum" "Embryo" "Heart" "Hypothalamus" "Ileum" "Jejunum" "Kidney" "Leukocytes" "Liver" "Lung" "Macrophage" "Muscle" "Ovary" "Oviduct" "Pituitary" "Retina" "Skin" "Small_intestine" "Spleen" "Testis" "Thymus" "Trachea")

mkdir -p ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 1. output all nominal SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results
rm -f -r ${dir_output}/nominal_combined_files.txt
nominal_combined_files=(`ls ${dir_nominal}/*/*.cis_qtl_pairs.*.parquet.csv.gz`)
for l in ${nominal_combined_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_combined_files.txt
}
done
python3 combine_signif_pairs_tjy.py ${dir_output}/nominal_combined_files.txt nominal_pairs -o ${dir_output}
#> output file: nominal_pairs.combined_signifpairs.txt.gz

### 2. extract all nominal pairs from nominal results for each tissue
for tis_i in {0..27}
do
{
    tissue=${tis_names[tis_i]}
    nominal_files2=(`ls ${dir_nominal}/${tissue}/ChickenGTEx.${tissue}.cis_qtl_pairs.*.parquet.csv.gz`)
    rm -f ${dir_output}/${tissue}.nominal_files2.txt
    for l in ${nominal_files2[*]}
    do
    {
        echo ${l} >> ${dir_output}/${tissue}.nominal_files2.txt
    }
    done
    # extract_pairs
    python3 extract_pairs_tjy.py ${dir_output}/${tissue}.nominal_files2.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${tissue}_nominal_pairs -o ${dir_output}
    #> output file: *_nominal_pairs.extracted_pairs.txt.gz
}
done
wait

### 3. prepare random SNP-gene pairs for MashR
nominal_pairs_files=(`ls ${dir_output}/*_nominal_pairs.extracted_pairs.txt.gz`)
rm -f ${dir_output}/nominal_pairs_files.txt
for l in ${nominal_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_pairs_files.txt
}
done
# MashR format file (z-score)
python3 mashr_prepare_input.py ${dir_output}/nominal_pairs_files.txt nominal_pairs.${subset_size}_subset -o ${dir_output} --only_zscore --dropna --subset $subset_size --seed 9823
zcat ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz | sed -e 's/.nominal_pairs_zval//g' | gzip > ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz
mv ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz
