library(data.table)
library(tidyverse)
library(RColorBrewer)
library(arrow, warn.conflicts = FALSE)


phenotype1="ENSGALG00000000249"
t="Liver"
c="Erythrocytes"
##### load snp map
map=fread(paste0("t,"/ChickenGTEx.", t,".bim"), nThread=48, tmpdir="./Temp") %>% as.data.frame
names(map)=c("CHR", "variant_id", "pos0", "POSITION", "Alt", "Ref")

######################### qtl1
#### load data
lead_df1 = fread(paste0("t,"/ChickenGTEx.", t,".cis_qtl.fdr.txt"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype1)
merge(lead_df1, map, by = "variant_id") -> lead_df1
chr = as.integer(lead_df1$CHR)
lead_snp1=as.character(lead_df1$variant_id)
top_snp="rs14303039"
qtl_df1=fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl_pairs.", chr,".parquet.csv.gz"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == phenotype1)
merge(qtl_df1, map, by = "variant_id") -> qtl_df1

#### determine signif threshold
threshold1=lead_df1[lead_df1$phenotype_id == phenotype1,"pval_nominal_threshold"] %>% as.numeric

#### compute LD 
if (!dir.exists("Temp")){
  dir.create("Temp")
}
ld_file=tempfile(pattern = "ld_file.", "Temp")
bfile=paste0("t,"/ChickenGTEx.", t)
#cmd=paste0("plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snp ", lead_snp1, " --ld-window-kb 5000 --ld-window 10000 --ld-window-r2 0 --r2 --out ", ld_file)
cmd=paste0("plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snps ", top_snp, " rs315662005 --ld-window-kb 10000  --ld-window 99999 --ld-window-r2 0 --r2 --out ", ld_file)
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
bfile=paste0("t,"/ChickenGTEx.", t)
cell = read_parquet(paste0("t,"/", c,"/", c,".cis_qtl_pairs.", chr,".parquet"),as_data_frame = TRUE) %>% filter(phenotype_id == phenotype1)
bim=fread(paste0("t,"/ChickenGTEx.", t,".bim"), header=F)
names(bim)=c("chr", "variant_id", "pos0", "pos", "alt", "red")
merge(cell, bim, by = "variant_id") -> cell
#top_snp=as.character(cell[which.min(cell$pval_gi),"variant_id"])
ld_file=tempfile(pattern = "ld_file.", "Temp")
cmd=paste0("plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snps ", top_snp, " rs315662005 --ld-window-kb 10000  --ld-window 99999 --ld-window-r2 0 --r2 --out ", ld_file)
system(cmd)
ld=fread(paste0(ld_file, ".ld"), header=T) %>% as.data.frame %>% select(SNP_B, R2)
names(ld)=c("SNP", "LD")
ld$LD=round(ld$LD, 1)
ld=rbind(ld, data.frame(SNP=top_snp, LD="1"))
merge(cell, ld, by.x="variant_id", by.y="SNP", all.x=T) ->  cell
cell$LD=ifelse(is.na(cell$LD), 0, cell$LD)
cell$LD=as.numeric(cell$LD)
ylim_value=ceiling(-log10(min(qtl_df1$pval_nominal, cell$pval_gi)))

p1=ggplot()+
      geom_point(qtl_df1 %>% filter(POSITION > 3650000 & POSITION < 3950000), mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal), color=LD), size=0.6, shape = 19, alpha=0.5)+
      geom_point(qtl_df1 %>% filter(variant_id == top_snp), mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal)), fill="red", color="red", shape=23, size=1.5)+
	  geom_point(qtl_df1 %>% filter(variant_id == "rs315662005"), mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal)), fill="red", color="red", shape=17, size=1.5)+
      geom_hline(yintercept=-log10(threshold1), linetype="dashed", color = "red")+
      #ylim(0, ylim_value)+
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
      geom_point(cell %>% filter(pos > 3650000 & pos < 3950000), mapping=aes(x = pos/1000000, y = -log10(pval_gi), color=LD), size=0.6, shape = 19, alpha=0.5)+
      geom_point(cell %>% filter(variant_id == top_snp), mapping=aes(x = pos/1000000, y = -log10(pval_gi)), fill="red", color="red", shape=23, size=1.5)+
	  geom_point(cell %>% filter(variant_id == "rs315662005"), mapping=aes(x = pos/1000000, y = -log10(pval_gi)), fill="red", color="red", shape=17, size=1.5)+
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

library(ggpubr)
p=ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend="right")
  ggsave(paste0("exampples.pdf"), p, width=4.5, height=3.3)