#!/bin/bash
source ~/.bashrc

dir_output="./"
dir_nominal="/group/zhougrp2/dguan/90_updates/02_eQTLs/04_results/"

tis_names=("Adipose" "Blood" "Brain" "Bursa" "Cecum" "Cerebellum" "Duodenum" "Embryo" "Heart" "Hypothalamus" "Ileum" "Jejunum" "Kidney" "Leukocytes" "Liver" "Lung" "Macrophage" "Muscle" "Ovary" "Oviduct" "Pituitary" "Retina" "Skin" "Small_intestine" "Spleen" "Thymus" "Testis" "Trachea")

mkdir -p  ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results
rm -f ${dir_output}/permutation_files.txt
perm_files=(`ls ${dir_nominal}/*/ChickenGTEx.*.cis_qtl.fdr.txt`)  # permutation results
for l in ${perm_files[*]}
do
{
    echo ${l} >> ${dir_output}/permutation_files.txt
}
done
python3 combine_signif_pairs_tjy.py ${dir_output}/permutation_files.txt strong_pairs -o ${dir_output}
#> output file: strong_pairs.combined_signifpairs.txt.gz

### 2. extract top pairs from nominal results for each tissue
for tis_i in {0..27}
do
{
    tissue=${tis_names[tis_i]}
    nominal_files=(`ls ${dir_nominal}/${tissue}/ChickenGTEx.${tissue}.cis_qtl_pairs.*.parquet.csv.gz`)
    rm -f ${dir_output}/${tissue}.nominal_files.txt
    for l in ${nominal_files[*]}
    do
    {
        echo ${l} >> ${dir_output}/${tissue}.nominal_files.txt
    }
    done
    # extract_pairs
    python3 extract_pairs_tjy.py ${dir_output}/${tissue}.nominal_files.txt ${dir_output}/strong_pairs.combined_signifpairs.txt.gz ${tissue} -o ${dir_output}
    #> output file: *.extracted_pairs.txt.gz
} &
done
wait

### 3. prepare strong SNP-gene pairs for MashR
strong_pairs_files=(`ls ${dir_output}/*.extracted_pairs.txt.gz`)
rm -f ${dir_output}/strong_pairs_files.txt
for l in ${strong_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/strong_pairs_files.txt
}
done
# MashR format file (z-score)
python3 mashr_prepare_input.py ${dir_output}/strong_pairs_files.txt strong_pairs -o ${dir_output} --only_zscore
zcat ${dir_output}/strong_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > ${dir_output}/strong_pairs.temp.txt.gz
mv ${dir_output}/strong_pairs.temp.txt.gz ${dir_output}/strong_pairs.MashR_input.txt.gz