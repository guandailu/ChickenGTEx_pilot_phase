library(tidyverse)
library(data.table)


color_df = fread("Tissues_for_eqtl_mapping.txt") %>% as.data.frame
names(color_df)=c("tissues", "samples", "colors")
tissues=color_df$tissues
mydf =data.frame()
for (t in tissues){
    df = fread(paste0(t,".cis_independent_qtl.txt.gz")) %>% as.data.frame %>% mutate(Tissue = t)
    mydf = rbind(mydf, df)
}
mydf$indep = ifelse(mydf$rank <= 3, mydf$rank, ">=4")
mydf$tss = abs(mydf$tss_distance)
mydf$maf = ifelse(mydf$af <= 0.5, mydf$af, 1-mydf$af)
mydf %>% select(indep, tss, maf) -> plot_df

ndf =data.frame()
for (t in tissues){
    df = fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl.txt")) %>% as.data.frame %>% mutate(Tissue = t) %>% filter(is_eGene == "FALSE")
    ndf = rbind(ndf, df)
}
nonegnenes=fread("non-eGenes.list", header=F) %>% pull(V1) %>% as.vector
ndf %>% subset(phenotype_id %in% nonegnenes) -> ndf
ndf$indep = 0
ndf$tss = abs(ndf$tss_distance)
ndf$maf = ifelse(ndf$af <= 0.5, ndf$af, 1-ndf$af)
ndf %>% select(indep, tss, maf) -> plot_df0

plot_df=rbind(plot_df, plot_df0)

p=ggplot(plot_df, aes(x=reorder(indep, tss/1000, FUN = median), y=tss/1000, fill=factor(indep, levels=c("0", "1", "2", "3", ">=4"))))+
    #geom_density(adjust=1.5, alpha=.4) +
    geom_boxplot(alpha=0.6)+
    xlab("Genes with # indep SNPs")+
    ylab("Distance to TSS (kb)")+
    #coord_flip()+
    #facet_wrap(vars(indep), nrow = 5, strip.position = "right")+
    scale_fill_manual(values=c("#808080", "#E9212C", "#01844F", "#7195C5", "#A64294"))+
    theme_classic(base_line_size = 0.4)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.text.x = element_text(color="black", size=8),
          axis.title.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=8),
          axis.title.y = element_text(color="black", size=10),
          #strip.background = element_rect(color="black", fill="#C0C0C0", size=0.2, linetype="solid"),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("TSS_distance_by_indep.boxplot.with_nonegene.pdf", p, width=3, height=3)
