### load required packages
load_packages <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)
  
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}

load_packages("data.table")
load_packages("tidyverse")

### parse parameters
ARGS <- commandArgs(trailingOnly = TRUE)
tspex = ARGS[1] # where tspex is install (tspex source: https://github.com/apcamargo/tspex)
expr_matrix = ARGS[2] # TPM expression matrix (row is gene id, coloumn is samples), colnames with sample ID is required (header: GeneID,sample1, sample2,...)
meta_tab = ARGS[3] # A metadata with tissue types (header: SampleID, Tissue)
prefix = ARGS[4] # output prefix

### Create temporary directory
temp_dir="Temp"
if (!dir.exists(temp_dir)){
  dir.create(temp_dir)
}

### Prepare tpm matrix
expr = as.data.frame(fread(expr_matrix, header=T))
meta = as.data.frame(fread(meta_tab, header=T))
rownames(expr)=as.vector(expr[,1])
expr = as.data.frame(t(expr))
expr$SampleID = rownames(expr)
expr= merge(meta, expr, by = "SampleID") 
gather(., GeneID, TPM, 3:ncol(expr), factor_key=TRUE) -> expr
expr %>% group_by(GeneID, Tissue) %>% summarise(mTPM=median(TPM)) %>% spread(Tissue, mTPM) %>% as.data.frame() -> expr
rownames(expr)=expr$GeneID
expr$GeneID=NULL
tpm_output=tempfile(pattern = paste0(prefix, "."), tmpdir = temp_dir)
fwrite(expr, tpm_output, sep="\t", row.names=T)


### run tspex
idx=c("counts", "tau", "gini", "simpson", "shannon_specificity", "roku_specificity", "tsi", "zscore", "spm", "spm_dpm", "js_specificity", "js_specificity_dpm")
for (x in idx){
 cmd=paste0(tspex, " --log ", tpm_output," ", prefix,".", x,".txt ", x)
 system(cmd)
}


### Remove temporary file
unlink(temp_dir, recursive = TRUE)