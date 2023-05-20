# October 21,2022
# Version 0
# Author: Dailu Guan
### load required packages
load_packages <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)
  
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}

### parse parameters
load_packages("argparser")
arguments <- arg_parser("****** The script is used for running fastenloc")
arguments <- add_argument(arguments, "mainfunc", default="LDblockCalc", help="Functions (inclduing LDblockCalc/CalcZscore/CalcPIP/DAP2enloc/fastenloc) to be run")
arguments <- add_argument(arguments, "--file", short="-f", help="prefix of plink file")
arguments <- add_argument(arguments, "--gwas", short="-g", help="GWAS file name")
arguments <- add_argument(arguments, "--snp", short="-id", help="SNP ID column name of the GWAS file")
arguments <- add_argument(arguments, "--is_zscore", short="-z", default='FALSE', help="z_score column name of the GWAS file")
arguments <- add_argument(arguments, "--beta", short="-b", help="beta column name of the GWAS file")
arguments <- add_argument(arguments, "--se", short="-s", help="standard error column name in GWAS file")
arguments <- add_argument(arguments, "--plink", short="-p", help="/path/to/plink")
arguments <- add_argument(arguments, "--torus", short="-r", help="/path/to/torus")
arguments <- add_argument(arguments, "--fastenloc", short="-e", help="/path/to/fastenloc")
arguments <- add_argument(arguments, "--dap2enloc", short="-d2e", help="/path/to/summarize_dap2enloc.pl, this script is provided by fastenloc package")
arguments <- add_argument(arguments, "--fastenloc_ann", short="-fa", help="fastenloc annotation file")
arguments <- add_argument(arguments, "--tissue", short="-ts", help="tissue where eQTLs are identified")
arguments <- add_argument(arguments, "--out_prefix", short="-o", help="output file prefix")
arguments <- add_argument(arguments, "--thread", short="-j", default='1',help="number of parallel threads for analysis")
arguments <- add_argument(arguments, "--tempdir", short="-t", default='Temp', help="/path/to/temporary/folder")
argv <- parse_args(arguments)

load_packages("data.table")
load_packages("tidyverse")

### Create temporary directory
if (!dir.exists(argv$tempdir)){
    dir.create(argv$tempdir)
}

### Create LD bolck by gwas
plink=argv$plink
if (file.exists(paste0(argv$file, ".bim"))){
	file_inputs=paste0("--bfile ", argv$file)
	bim=paste0(argv$file, ".bim")
}else{
	file_inputs=paste0("--file ", argv$file)
	bim=paste0(argv$file, ".map")
}
bim_df=fread(bim, header=F) %>% as.data.frame()
if (argv$mainfunc == "LDblockCalc"){
    cat("Running LDblockCalc...\n\n")
    ld_prefix=str_replace(argv$file,".*/","")
    cmd=paste0(plink, " ", file_inputs, " --chr-set 39 --blocks no-pheno-req --blocks-max-kb 1000 --make-founders --out ", ld_prefix)
    system(cmd)
}

if (argv$mainfunc == "CalcZscore"){
    gwas_zscore_file=paste0(argv$out_prefix, ".zscore.txt.gz")
    if (!file.exists(gwas_zscore_file)){
	bfile_prefix=str_replace(argv$file, ".*/","")
        blocks=fread(paste0(bfile_prefix, ".blocks.det"), header=T) %>% as.data.frame()
	blocks$line_num=1:nrow(blocks)
        data.frame(line_num=rep(blocks$line_num, blocks$NSNPS), SNP=unlist(strsplit(blocks$SNPS, "\\|"))) -> mydf
    	bim = subset(bim_df, ! V2 %in% as.vector(unique(mydf$SNP)) )
    	rbind(mydf,data.frame(line_num=(max(mydf$line_num)+1):((max(mydf$line_num)+1)+nrow(bim)-1), SNP=bim$V2) ) -> mydf
    	gwas=fread(argv$gwas, header=T) %>% as.data.frame()
    	if (argv$is_zscore == "TRUE"){
            zscore=gwas[,"z"]
    	}else{
	        beta=argv$beta
	        se=argv$se
	       zscore=gwas[,beta]/gwas[,se]
	    }
	    snp=argv$snp
	    data.frame(SNP=gwas[,snp],zscore=zscore) -> gwas
	    merge(mydf,gwas, by="SNP", all=T) -> mydf
	    mydf = na.omit(mydf)
	    mydf = select(mydf, SNP, line_num, zscore)
	    fwrite(mydf, gwas_zscore_file, row.names=F, sep="\t", col.names=F, compress="gzip")
	}
}

### Calculate posterior inclusion probability
if (argv$mainfunc == "CalcPIP"){
    cat("Running PIP calculation...\n\n")
    gwas_zscore_file=paste0(argv$out_prefix, ".zscore.txt.gz")
    if (!file.exists(paste0(argv$gwas, ".pip.txt.gz"))){
        cmd=paste0(argv$torus, " -d ", gwas_zscore_file, " --load_zval -dump_pip ", argv$gwas, ".pip.txt")
        system(cmd)
        cmd=paste0("gzip -f ",argv$gwas, ".pip.txt")
        system(cmd)
    }
}

### DAP to fastenloc
if (argv$mainfunc == "DAP2enloc"){
    ### This was done by honghao zhong
    summarize_dap2enloc=argv$dap2enloc
    cmd=paste0(summarize_dap2enloc, " -dir ", dap_rst_dir, " -vcf ", snp_vcf_file, "-tissue", argv$tissue, " | gzip - > ", argv$tissue, ".fastenloc.eqtl.annotation.vcf.gz")
}

### Run fastenloc
if (argv$mainfunc == "fastenloc"){
    cat("Running fastenloc...\n\n")
    total_variants=fread(paste0(argv$gwas, ".pip.txt.gz"), header=F)
    length(unique(total_variants[,1])) -> total_variants
    #total_variants=nrow(bim_df)
    tissue=argv$tissue
    cmd=paste0(argv$fastenloc, " -eqtl ", argv$fastenloc_ann, " -gwas ", argv$gwas, ".pip.txt.gz -total_variants ", total_variants, " -thread ", argv$thread, " -prefix ", argv$out_prefix, ".", argv$tissue,".fastenloc.txt")
    system(cmd)
}
