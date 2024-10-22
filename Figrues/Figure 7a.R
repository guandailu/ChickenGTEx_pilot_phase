library(data.table)
library(tidyverse)

readRDS("chicken_TWAS_BW.rds") -> df
df %>% group_by(chr) %>% summarise(chr_len=max(tss_pos)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(df, ., by=c("chr"="chr")) %>% arrange(chr, tss_pos) %>% mutate(BPcum=tss_pos+tot) -> df
fdr_threshold=df %>% filter(fdr_ref > 0.05) %>% slice(which.min(pvalue)) %>% pull(pvalue)
df_center = df %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
library(RColorBrewer)
col_vector = rep(c("#000000", "#A9A9A9"), 16)
library(ggrepel)
df$color_chr = as.character(paste0("chr", df$chr))
p=ggplot(df)+
    geom_point(mapping=aes(x = BPcum, y = -log10(pvalue), color=color_chr), size=1.2, shape = 19)+
    geom_point(df %>% subset(fdr_ref< 0.05), mapping=aes(x =BPcum, y = -log10(pvalue)), shape=18, color="#000080", fill="#000080", size=3)+
    geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color = "red")+    scale_color_manual(values=col_vector)+
    #scale_shape_manual(values=c(shapes))+
    #coord_flip()+
    xlab("")+
    ylab(expression("-log"[10]*italic(P)))+
    scale_x_continuous(breaks = c(df_center$center)[c(1:5, 7, 10, 15, 20, 31)], label = c(df_center$chr)[c(1:5, 7, 10, 15, 20, 31)])+
    theme_classic(base_size = 15, base_line_size = 0.6)+
    theme(legend.title=element_blank(), legend.position="none",legend.text = element_text(size=10),legend.key=element_blank(),legend.box.background = element_blank(),legend.background = element_blank(),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(2, 1, 2, 2))+
    guides(color=guide_legend(ncol=1)) +
    geom_text_repel(
        data = subset(df, fdr_ref < 0.05),
        aes(x= BPcum, y = -log10(pvalue), label =gene_name),
        direction= "x", size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
    )
ggsave("chicken_body_weight_within_species_TWAS.pdf", p, width=8, height=4)
library(qqman)
pdf("chicken_body_weight_within_species_TWAS.qqplot.pdf", width=3.5, height=4)
qq(df$pvalue)
dev.off()

readRDS("chicken_TWAS_metaTWAS.rds") -> df
p=ggplot(df)+
    geom_point(mapping=aes(x = BPcum, y = log10(pvalue_meta), color=color_chr), size=1.2, shape = 19)+
    geom_point(df %>% subset(fdr_meta< 0.05 & is_signif == "both_sig"), mapping=aes(x =BPcum, y =log10(pvalue_meta)), shape=18, color="#000080", fill="#000080", size=3)+
    geom_point(df %>% subset(fdr_meta< 0.05 & is_signif == "meta_novel_sig"), mapping=aes(x =BPcum, y =log10(pvalue_meta)), shape=1, color="red", size=2.5)+
    geom_point(df %>% subset(fdr_meta< 0.05 & is_signif == "meta_novel_sig" & gene %in% validated_genes), mapping=aes(x =BPcum, y =log10(pvalue_meta)), shape=21, color="red", fill="red", size=2.5)+
    geom_hline(yintercept=log10(fdr_threshold), linetype="dashed", color = "red")+    scale_color_manual(values=col_vector)+
    #scale_shape_manual(values=c(shapes))+
    #coord_flip()+
    xlab("")+
    ylab(expression("-log"[10]*italic(P)))+
    scale_x_continuous(position = "top", breaks = c(df_center$center)[c(1:5, 7, 10, 15, 20, 31)], label = c(df_center$chr)[c(1:5, 7, 10, 15, 20, 31)])+
    scale_y_break(c(-40, -70))+
    theme_classic(base_size = 15, base_line_size = 0.6)+
    theme(legend.title=element_blank(), legend.position="none",legend.text = element_text(size=10),legend.key=element_blank(),legend.box.background = element_blank(),legend.background = element_blank(),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(2, 1, 2, 2))+
    guides(color=guide_legend(ncol=1)) +
    geom_text_repel(
        data = subset(df, fdr_meta< 0.05 & is_signif == "meta_novel_sig" & gene %in% validated_genes),
        aes(x= BPcum, y = log10(pvalue_meta), label =gene_name),
        nudge_x = .15,
        box.padding = 0.5,
        nudge_y = 1,
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 30
    )
ggsave("chicken_body_weight_metaTWAS.pdf", p, width=8, height=4)
pdf("chicken_body_weight_metaTWAS.qqplot.pdf", width=3.5, height=4)
qq(df$pvalue_meta)
dev.off()