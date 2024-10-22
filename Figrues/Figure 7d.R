library(data.table)
library(tidyverse)

read.table("validated_genes.list", header=F) %>% pull(V1) -> validated_genes
readRDS("chicken_TWAS_BW.rds") -> metatwas

for (g in validated_genes){
gene_name=metatwas[metatwas$gene == g,"gene_name"]
phewas_file=paste0(gene_name, ".phewas.csv")
if (file.exists(phewas_file)){
df = fread(phewas_file)
df$fdr = p.adjust(df$`P-value`, method = "fdr")
df %>% group_by(Domain) %>% mutate(x_num=row_number()) %>% as.data.frame -> df
df %>% group_by(Domain) %>% summarise(chr_len=max(x_num)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(df, ., by=c("Domain"="Domain")) %>% arrange(Domain, x_num) %>% mutate(BPcum=x_num+tot) -> df
fdr_threshold=df %>% filter(fdr > 0.05) %>% slice(which.min(`P-value`)) %>% pull(`P-value`)
df_center = df %>% group_by(Domain) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
library(RColorBrewer)
n <- 28
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df$trait_lab = str_remove(df$Trait, "Impedance measures - ")
df$Domain = ifelse(df$Domain == "Environmental", "Environment", df$Domain)
library(ggrepel)
p=ggplot(df)+
    geom_point(mapping=aes(x = BPcum, y = -log10(`P-value`), color=Domain), size=0.9, shape = 19)+
    geom_point(df %>% subset(fdr< 0.05), mapping=aes(x =BPcum, y = -log10(`P-value`), color=Domain), shape=13, size=2)+
    geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color = "red")+    scale_color_manual(values=col_vector)+
    #scale_shape_manual(values=c(shapes))+
    #coord_flip()+
    xlab("")+
    ylab(expression("-log"[10]*italic(P)))+
    scale_x_continuous(breaks = c(df_center$center), label = NULL)+
    theme_classic(base_size = 15, base_line_size = 0.6)+
    theme(legend.title=element_blank(), legend.position="right",legend.text = element_text(size=10),legend.key=element_blank(),legend.box.background = element_blank(),legend.background = element_blank(),
          axis.text.x = element_text(color="black", size=10, angle=90, hjust=1, vjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(2, 1, 2, 2))+
    guides(color=guide_legend(ncol=1)) #+
    #geom_text_repel(
    #    data = subset(df, fdr < 0.0001),
    #    aes(x= BPcum, y = -log10(`P-value`), label =trait_lab),
    #    direction= "x", size = 5,
    #    box.padding = unit(0.35, "lines"),
    #    point.padding = unit(0.3, "lines")
    #)
ggsave(paste0(gene_name, ".phewas.pdf"), p, width=8, height=2)
}
}