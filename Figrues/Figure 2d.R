library(tidyverse)
library(data.table)


df1 = fread("/group/zhougrp2/dguan/90_updates/02_eQTLs/04_results/expression_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "eQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)
df2 = fread("/group/zhougrp2/dguan/90_updates/03_sQTLs/04_results/splicing_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "sQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)
df3= fread("/group/zhougrp2/dguan/90_updates/04_lncQTLs/04_results/lncRNA_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "lncQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)
df4= fread("/group/zhougrp2/dguan/90_updates/04_lncQTLs/04_results/lncRNA_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "lncQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)
df4= fread("/group/zhougrp2/dguan/90_updates/05_exQTLs/04_results/exon_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "exQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)
df5 = fread("/group/zhougrp2/dguan/90_updates/06_3aQTLs/04_results/APA_qtl_summary.txt") %>% as.data.frame %>% mutate(QTL = "3aQTL") %>% select(tissues, samples, detected_genes, egenes, QTL)


df = rbind(df1, df2, df3, df4, df5)
color_df = fread("/group/zhougrp2/dguan/90_updates/Tissues_with_avail_samples_for_eqtl_mapping.txt") %>% as.data.frame
names(color_df)=c("tissues", "samples", "colors")
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df$tissues)


p=ggplot(df)+
    geom_point(aes(x=samples, y= egenes/detected_genes, color=tissues), shape=19, size=3)+
    #geom_smooth(aes(x=samples, y= egenes/detected_genes), color="red", method="lm", se=F)+
    facet_wrap(vars(factor(QTL, levels=c("eQTL", "sQTL", "lncQTL", "exQTL", "3aQTL"))), nrow = 1)+
    ylim(0, 1)+
    xlab(expression("Sample size"))+
    ylab("#ePhenotype/#tested phenotypes")+
    scale_color_manual(values=colors)+
    scale_shape_manual(values=c(19, 9, 17, 8, 15))+
    theme_classic(base_line_size = 0.4)+
    theme(legend.title = element_blank(), legend.position="right",
          axis.text.x = element_text(color="black", size=8),
          axis.title.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=8),
          axis.title.y = element_text(color="black", size=10),
          strip.background = element_rect(color="black", fill="#C0C0C0", size=0.2, linetype="solid"),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+guides(color=guide_legend(nrow=3))
ggsave("sample_size_vs_xgenes.no-smooth.pdf", p, width=7.5, height=2.3)

