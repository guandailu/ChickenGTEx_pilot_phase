library(tidyverse)
library(data.table)

df1 = fread("Gallus_gallus.GRCg6a.102.exons.phastCons77way.tab", header=F) %>% as.data.frame
df1 %>% separate(., V1, into=c("Gene", "exon", "strand"), sep=":") -> df1
df1 %>% group_by(Gene) %>% summarize(phastcons=mean(V6)) -> df1
merge(mydf, df1, by.x="phenotype_id", by.y="Gene") %>% select(indep, phastcons) -> plot_df

nonegnenes=fread("non-eGenes.list", header=F) %>% as.data.frame
names(nonegnenes)="Gene"
merge(nonegnenes, df1, by="Gene") %>% mutate(indep=0) %>% select(indep, phastcons) -> nonegnenes
rbind(plot_df, nonegnenes) -> plot_df

p=ggplot(plot_df, aes(x=factor(indep, levels=c("1", "2", "3", ">=4", "0")), y=phastcons, fill=factor(indep, levels=c("0", "1", "2", "3", ">=4"))))+
    #geom_density(adjust=1.5, alpha=.4) +
    geom_boxplot(alpha=0.6)+
    xlab("eGene with # indep SNPs")+
    ylab("phastCons")+
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
ggsave("phastcons_by_indep.boxplot.with_nonegene.pdf", p, width=3, height=3)
