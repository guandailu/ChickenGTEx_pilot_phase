library(data.table)
library(tidyverse)
tissues=fread("Tissues_for_eqtl_mapping.txt", header=F) %>% pull(V1) %>% as.vector
mydf = data.frame()
for (t in tissues){
df = fread(paste0(t,"/ChickenGTEx.",t,".cis_qtl.txt")) %>% as.data.frame %>% select(phenotype_id, num_var) %>% mutate(Tissue = t)
mydf = rbind(mydf, df)
}
mydf = mydf %>% group_by(phenotype_id) %>% summarize(num_snps=as.integer(mean(num_var)))
p=ggplot(mydf)+
    geom_histogram(aes(x=num_snps), binwidth = 500) +
    geom_vline(xintercept = median(mydf$num_snps), linetype="dashed", color = "red", linewidth=1)+
    xlab("# SNPs tested/gene")+
    ylab("Number of genes")+
    ggplot2::theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=15),
          axis.text.y =element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("Accuracy/numSNPs_per_tested_gene.pdf", p, width=4, height=3.5) 
