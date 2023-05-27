#!/bin/bash

set -e

source ~/.bashrc

wdr=$1

tissues=("Adipose" "Blood" "Brain" "Bursa" "Cecum" "Cerebellum" "Duodenum" "Embryo" "Heart" "Hypothalamus" "Ileum" "Jejunum" "Kidney" "Leukocytes" "Liver" "Lung" "Macrophage" "Muscle" "Ovary" "Oviduct" "Pituitary" "Retina" "Skin" "Small_intestine" "Spleen" "Testis" "Thymus" "Trachea")

rm -f ${wdr}/signifpair_list_files.list
for t in ${tissues[@]};
do
    echo -e "${wdr}/../04_results/${t}/ChickenGTEx.${t}.cis_independent_qtl.txt.gz" >> ${wdr}/signifpair_list_files.list
done

python3 ${wdr}/step1.combine_signif_pairs.py \
    ${wdr}/signifpair_list_files.list \
	${wdr}/ChickenGTEx.SNP_map.txt \
	ChickenGTEx \
	--output_dir ${wdr}

for t in ${tissues[@]};
do 
    rm -f ${t}.nominal_files.txt
    for chr in {1..28} {30..33};
    do 
      echo ${wdr}/../04_results/${t}/ChickenGTEx.${t}.cis_qtl_pairs.${chr}.parquet.csv.gz >> ${t}.nominal_files.txt
    done
    python3 ${wdr}/step2.extract_pairs.py \
        ${t}.nominal_files.txt \
        ${wdr}/ChickenGTEx.combined_signifpairs.txt.gz \
        ChickenGTEx.${t} \
        --output_dir ${wdr}
done


rm -f ${wdr}/variant_gene_pair_files.list
for t in ${tissues[@]};
do 
    echo -e "${wdr}/ChickenGTEx.${t}.extracted_pairs.txt.gz" >> ${wdr}/variant_gene_pair_files.list
done

python3 ${wdr}/step3.metasoft_prepare_input.py \
    ${wdr}/variant_gene_pair_files.list \
	ChickenGTEx \
	--output_dir ${wdr} \
	--write_full


python3 ${wdr}/step4.run_metasoft.py \
    /home/dguan/bin/metasoft/Metasoft.jar \
	ChickenGTEx.metasoft_input.txt.gz \
	ChickenGTEx \
	--pvalue_table /home/dguan/bin/metasoft/HanEskinPvalueTable.txt \
	--seed 666666 \
	--output_dir ${wdr}

rm -rf ${wdr}/metasoft_output.list
echo -e "${wdr}/ChickenGTEx.metasoft.txt.gz" >> ${wdr}/metasoft_output.list
python3 ${wdr}/step5.metasoft_postprocess.py \
        ${wdr}/metasoft_output.list \
	${wdr}/variant_gene_pair_files.list \
	ChickenGTEx.postprocess \
	--output_dir ${wdr}