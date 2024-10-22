library(data.table)
library(tidyverse)
df = read.table("validated_genes_results.txt", header=T)
df %>% filter(pvalue < 0.05) %>% group_by(gene_name) %>%  slice(which.min(pvalue)) -> df
p= ggplot(df)+
      geom_bar(mapping=aes(x = reorder(gene_name, log10(pvalue)), y = -log10(pvalue)), stat="identity", fill="#92C5DE")+
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+
      xlab("")+
      ylab(expression("-log"[10]*italic(P)))+
      coord_flip()+
      theme_classic(base_size = 15, base_line_size = 0.6)+
      theme(legend.title=element_blank(), legend.position="none",legend.text = element_text(size=10),legend.key=element_blank(),legend.box.background = element_blank(),legend.background = element_blank(),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(2, 1, 2, 2))
ggsave("metaTWAS.validated_genes.pdf", p, width=3, height=3)