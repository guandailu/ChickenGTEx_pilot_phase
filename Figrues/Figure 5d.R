library(data.table)
library(tidyverse)
library(RColorBrewer)

ARGS <- commandArgs(trailingOnly = TRUE)
variant_id=ARGS[1] #e.g. "rs741492341"
gene=ARGS[2] #"e.g. ENSGALG00000026471"
t=ARGS[3]  #"e.g. Jejunum"

##### load snp map
map=fread(paste0("t,"/ChickenGTEx.", t,".bim"), nThread=48, tmpdir="./Temp") %>% as.data.frame
names(map)=c("CHR", "variant_id", "pos0", "POSITION", "Alt", "Ref")

#### load data
chr = as.integer(map[map$variant_id == variant_id,"CHR"])
lead_df = fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl.txt"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == gene)
merge(lead_df, map, by = "variant_id") -> lead_df
eqtl=fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl_pairs.", chr,".parquet.csv.gz"), nThread=48) %>% as.data.frame %>% filter(phenotype_id == gene)
merge(eqtl, map, by = "variant_id") -> eqtl
threshold=lead_df[lead_df$phenotype_id == gene,"pval_nominal_threshold"] %>% as.numeric

#### region
left_boundary = as.integer(min(eqtl$POSITION))
right_boundary = as.integer(max(eqtl$POSITION))

#### compute LD 
if (!dir.exists("Temp")){
  dir.create("Temp")
}
ld_file=tempfile(pattern = "ld_file.", "Temp")
bfile=paste0("t,"/ChickenGTEx.", t)
cmd=paste0(plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --ld-snp ", variant_id, " --ld-window-kb 5000 --ld-window 10000 --ld-window-r2 0 --r2 --out ", ld_file)
system(cmd)
ld=fread(paste0(ld_file, ".ld"), header=T) %>% as.data.frame %>% select(SNP_B, R2)
names(ld)=c("SNP", "LD")
ld=rbind(ld, data.frame(SNP=variant_id, LD="1"))
merge(eqtl, ld, by.x="variant_id", by.y="SNP", all.x=T) ->  eqtl
eqtl$LD=ifelse(is.na(eqtl$LD), 0, eqtl$LD)
eqtl$LD=as.numeric(eqtl$LD)

p4=ggplot()+
      geom_point(eqtl, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal), color=LD), size=0.9, shape = 19)+
      geom_point(lead_df, mapping=aes(x = POSITION/1000000, y = -log10(pval_nominal)), color="red", fill="red", shape=23, size=2)+
      geom_hline(yintercept=-log10(threshold), linetype="dashed", color = "red")+
      #scale_color_manual(values=c(ld_colors))+
      xlab("Position (Mb)")+
      ylab(expression("-log"[10]*italic(P)))+
      #scale_y_continuous(breaks=seq(0,20,10), limits = c(0, 22))+
      scale_color_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")), limits = c(0, 1), breaks = c(seq(0, 1, 0.2)), label = c(seq(0, 1, 0.2))) + 
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.title=element_blank(),
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(1, 1, 1, 1))
ggsave(paste0(t, ".", gene, ".", variant_id, ".pdf"), p4, width=5, height=3)
