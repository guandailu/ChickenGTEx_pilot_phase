library(data.table)
library(tidyverse)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue = ARGS[1]
cat("Processing tissue: ", tissue, "...\n")
signif_pair <- function(tissue){
  fdrdf=fread(paste0(tissue,"/ChickenGTEx.", tissue,".cis_qtl.fdr.txt"))
  fdrdf=fdrdf[,c("phenotype_id", "pval_nominal_threshold")]
  for (chr in c(1:28, 30:33)){
    df=fread(paste0(tissue, "/ChickenGTEx.",tissue,".cis_qtl_pairs.", chr,".parquet.csv.gz"))
    DF = merge(df, fdrdf, by="phenotype_id")
    DF=filter(DF, pval_nominal < pval_nominal_threshold)
    if (nrow(DF) > 0){
      fwrite(DF, paste0(tissue,"/ChickenGTEx.", tissue,".signifpairs",chr,".txt.gz"), sep = "\t", row.names=F, compress="gzip", quote=F)
    }else{
      cat("Chromosome ", chr, "does not have significant QTLs...\n")
    }
    rm(df)
    rm(DF)
  }
}
signif_pair(tissue)
DF=data.frame()
for (chr in c(1:28, 30:33)){
  signif_file=paste0(tissue,"/ChickenGTEx.", tissue,".signifpairs",chr,".txt.gz")
  if (file.exists(signif_file)){
    df = fread(signif_file, sep = "\t", header=T)
    DF=rbind(DF, df)
    system(paste0("rm -f ", signif_file))
  }
}
fwrite(DF, paste0(tissue,"/ChickenGTEx.", tissue,".signifpairs.txt.gz"), sep = "\t", row.names=F, compress="gzip", quote=F)

cat("Done!\n")
