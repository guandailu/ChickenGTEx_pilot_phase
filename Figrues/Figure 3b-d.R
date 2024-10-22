library(data.table)
library(tidyverse)
library(RColorBrewer)
ARGS <- commandArgs(trailingOnly = TRUE)
n=as.integer(ARGS[1])

load("molQTL_coloc.eqtl_apaqtl_coloc_results.Rdata")
#n=1
phenotype1=plot_df[n, "phenotype_id1"]
phenotype2=plot_df[n, "phenotype_id2"]
t=plot_df[n, "Tissue"]

##### load snp map
map=fread(paste0("/group/zhougrp2/dguan/90_updates/02_eQTLs/01_genotypes/", t,"/ChickenGTEx.", t,".bim"), nThread=48, tmpdir="./Temp") %>% as.data.frame
names(map)=c("CHR", "variant_id", "pos0", "POSITION", "Alt", "Ref")


######################### qtl1
#### load data
lead_df1 = fread(paste0("/group/zhougrp2/dguan/90_updates/02_eQTLs/04_results/", t,"/ChickenGTEx.", t,".cis_qtl.fdr.txt"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype1)
merge(lead_df1, map, by = "variant_id") -> lead_df1
chr = as.integer(lead_df1$CHR)
lead_snp1=as.character(lead_df1$variant_id)
qtl_df1=fread(paste0("/group/zhougrp2/dguan/90_updates/02_eQTLs/04_results/", t,"/ChickenGTEx.", t,".cis_qtl_pairs.", chr,".parquet.csv.gz"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype1)
merge(qtl_df1, map, by = "variant_id") -> qtl_df1

#### determine signif threshold
threshold1=lead_df1[lead_df1$phenotype_id == phenotype1,"pval_nominal_threshold"] %>% as.numeric

#### compute LD 
if (!dir.exists("Temp")){
  dir.create("Temp")
}
ld_file=tempfile(pattern = "ld_file.", "Temp")
bfile=paste0("/group/zhougrp2/dguan/90_updates/02_eQTLs/01_genotypes/", t,"/ChickenGTEx.", t)
cmd=paste0("/share/apps/plink-1.90/plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snp ", lead_snp1, " --ld-window-kb 5000 --ld-window 10000 --ld-window-r2 0 --r2 --out ", ld_file)
system(cmd)
ld=fread(paste0(ld_file, ".ld"), header=T) %>% as.data.frame %>% select(SNP_B, R2)
names(ld)=c("SNP", "LD")
#ld=rbind(ld, data.frame(SNP=variant_id, LD="1"))
ld$LD=as.numeric(ld$LD)
merge(qtl_df1, ld, by.x="variant_id", by.y="SNP", all.x=T) ->  qtl_df1
qtl_df1$LD=ifelse(is.na(qtl_df1$LD), 0, qtl_df1$LD)
system(paste0("rm -rf ", ld_file, ".ld"))


######################### qtl2
#### load data
lead_df2 = fread(paste0("/group/zhougrp2/dguan/90_updates/06_3aQTLs/04_results/", t,"/ChickenGTEx.", t,".cis_qtl.fdr.txt"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype2)
merge(lead_df2, map, by = "variant_id") -> lead_df2
chr = as.integer(lead_df2$CHR)
lead_snp2=as.character(lead_df2$variant_id)
qtl_df2=fread(paste0("/group/zhougrp2/dguan/90_updates/06_3aQTLs/04_results/", t,"/ChickenGTEx.", t,".cis_qtl_pairs.", chr,".parquet.csv.gz"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype2)
merge(qtl_df2, map, by = "variant_id") -> qtl_df2

#### determine signif threshold
threshold2=lead_df2[lead_df2$phenotype_id == phenotype2,"pval_nominal_threshold"] %>% as.numeric

#### compute LD 
if (!dir.exists("Temp")){
  dir.create("Temp")
}
ld_file=tempfile(pattern = "ld_file.", "Temp")
bfile=paste0("/group/zhougrp2/dguan/90_updates/06_3aQTLs/01_genotypes/", t,"/ChickenGTEx.", t)
cmd=paste0("/share/apps/plink-1.90/plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snp ", lead_snp2, " --ld-window-kb 5000 --ld-window 10000 --ld-window-r2 0 --r2 --out ", ld_file)
system(cmd)
ld=fread(paste0(ld_file, ".ld"), header=T) %>% as.data.frame %>% select(SNP_B, R2)
names(ld)=c("SNP", "LD")
#ld=rbind(ld, data.frame(SNP=variant_id, LD="1"))
ld$LD=as.numeric(ld$LD)
merge(qtl_df2, ld, by.x="variant_id", by.y="SNP", all.x=T) ->  qtl_df2
qtl_df2$LD=ifelse(is.na(qtl_df2$LD), 0, qtl_df2$LD)
system(paste0("rm -rf ", ld_file, ".ld"))

ylim_value=ceiling(-log10(min(qtl_df1$pval_nominal, qtl_df2$pval_nominal)))
p1=ggplot()+
      geom_point(qtl_df1, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal), color=LD), size=0.9, shape = 19)+
      geom_point(lead_df1, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal)), color="red", shape=18, size=4)+
      geom_hline(yintercept=-log10(threshold1), linetype="dashed", color = "red")+
      ylim(0, ylim_value)+
      #scale_color_manual(values=c(ld_colors))+
      xlab("")+
      ylab(expression("-log"[10]*italic(P)))+
      #scale_y_continuous(breaks=seq(0,20,10), limits = c(0, 22))+
      scale_color_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")), limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) + 
      theme_classic(base_size = 14, base_line_size = 0.5)+
      theme(legend.title=element_blank(),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(1, 1, 1, 1))


p2=ggplot()+
      geom_point(qtl_df2, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal), color=LD), size=0.9, shape = 19)+
      geom_point(lead_df2, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal)), color="red", shape=18, size=4)+
      geom_hline(yintercept=-log10(threshold2), linetype="dashed", color = "red")+
      ylim(0, ylim_value)+
      #scale_color_manual(values=c(ld_colors))+
      xlab("")+
      ylab(expression("-log"[10]*italic(P)))+
      #scale_y_continuous(breaks=seq(0,20,10), limits = c(0, 22))+
      scale_color_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")), limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) + 
      theme_classic(base_size = 14, base_line_size = 0.5)+
      theme(legend.title=element_blank(),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(1, 1, 1, 1))

if (!dir.exists("molQTL_coloc_example_plots")){
  dir.create("molQTL_coloc_example_plots")
}

library(ggpubr)
p=ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend="right")
  ggsave(paste0("molQTL_coloc_example_plots/coloc.eqtl_aoaqtl.",t, ".", phenotype1, ".", phenotype2,".pdf"), p, width=4, height=5)