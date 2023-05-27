library(data.table)
library(tidyverse)
ARGS <- commandArgs(trailingOnly = TRUE)
mol_pheno=ARGS[1] # eQTL sQTL lncQTL exQTL 3aQTL
gwas_file_prefix=ARGS[2] #  gwas file name

output_dir=paste0(mol_pheno, "_results")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

if (mol_pheno == "eQTL"){
    pheno="mRNA"
}else if (mol_pheno == "sQTL"){
    pheno="splicing"
}else if (mol_pheno == "lncQTL"){
    pheno="lncRNA"
}else if (mol_pheno == "exQTL"){
    pheno="exon"
}else if (mol_pheno == "3aQTL"){
    pheno="APA"
}

#tissues=c("Adipose", "Blood", "Brain", "Bursa", "Cecum", "Cerebellum", "Duodenum", "Embryo", "Heart", "Hypothalamus", "Ileum", "Jejunum", "Kidney", "Leukocytes", "Liver", "Lung", "Macrophage", "Muscle", "Ovary", "Oviduct", "Pituitary", "Retina", "Skin", "Small_intestine", "Spleen", "Testis", "Thymus", "Trachea")
tissues=c("Small-intestine")
if (pheno=="APA"){
    tissues=tissues[! tissues %in% c("Cerebellum", "Duodenum", "Spleen")]
}else{
    tissues=tissues
}

for (i in 1:length(tissues)){
  tissue=tissues[i]
  output_file=paste0(mol_pheno,"_results/", tissue, ".", gwas_file_prefix, ".twas.txt")
  if (file.exists(output_file)){
     file.remove(output_file)
  }
  cmd1=paste0("python3 /home/dguan/bin/MetaXcan/software/SPrediXcan.py --model_db_path /group/zhougrp2/dguan/41_GWAS_integrated/Hyline_GWAS/TWAS/db/",mol_pheno,"/ChickenGTEx_V1_", tissue, "_",pheno,"_ElasticNet_model_filtered_signif.db --covariance /group/zhougrp2/dguan/41_GWAS_integrated/Hyline_GWAS/TWAS/db/",mol_pheno,"/ChickenGTEx_V1_", tissue, "_",pheno,"_ElasticNet_models_covariances.gz --gwas_folder /group/zhougrp2/dguan/41_GWAS_integrated/Hyline_GWAS/Selected_GWAS --gwas_file_pattern ", gwas_file_prefix, "$ --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column b --pvalue_column p --output_file ", output_file)
  system(cmd1)
  
  cmd2=paste0("Rscript TWAS_plot.R ",mol_pheno, " ",mol_pheno,"_results/", tissue, ".", gwas_file_prefix, ".twas.txt ", mol_pheno,"_results/", tissue, ".", gwas_file_prefix, ".twas")
  system(cmd2)
}
