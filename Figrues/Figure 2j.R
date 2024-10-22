suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

df = fread("gene_ortholog.txt", header=T) %>% as.data.frame
tissues=fread("Tissues_for_eqtl_mapping.txt", header=F) %>% pull(V1)

egene_df=data.frame()
for (tiss in tissues){
  input_file=paste0(tiss,"/ChickenGTEx.", tiss,".cis_qtl.txt")
  tdf=fread(input_file, header=T) %>% as.data.frame %>% mutate(tissue=tiss) %>% select(phenotype_id, is_eGene)
  egene_df=rbind(egene_df, tdf)
}
egene_df$value = ifelse(egene_df$is_eGene == "TRUE", 1, 0)
egene_df %>% group_by(phenotype_id) %>% summarize(is_eGene = sum(value)) -> egene_df
egene_df$value = ifelse(egene_df$is_eGene == 0, "non-eGene", "eGene")
merge(egene_df, df, by.x="phenotype_id", by.y="Gene stable ID") -> df
df[df==""]<-NA
df[is.na(df)] <- "chicken_specific"
df$Type = ifelse(df$`Cow homology type` == "chicken_specific" & df$`Pig homology type` == "chicken_specific" & df$`Human homology type` == "chicken_specific", "Non-orthologs", ifelse(df$`Cow homology type` == "ortholog_one2one" & df$`Pig homology type` == "ortholog_one2one" & df$`Human homology type` == "ortholog_one2one", "Orthologs", "Other"))

df1 = fread("Gallus_gallus.GRCg6a.102.exons.phastCons77way.tab", header=F) %>% as.data.frame
df1 %>% separate(., V1, into=c("Gene", "exon", "strand"), sep=":") -> df1
df1 %>% group_by(Gene) %>% summarize(phastcons=mean(V6)) -> df1
merge(df1, df, by.x="Gene", by.y="phenotype_id") -> df

df %>% filter(Type != "Other") %>% select(Gene, Type, value, phastcons) -> df
#df %>% gather(type, ortholog, `Cow homology type`, `Human homology type`, `Pig homology type`) -> df
#df$ortholog = str_replace(df$ortholog, "ortholog_","")
#df$ortholog = factor(df$ortholog, levels=c("chicken_specific", "one2one", "one2many", "many2many"))
#df$type = str_replace(df$type, "Cow homology type","chicken vs cattle") %>% str_replace(., "Pig homology type","chicken vs pig") %>% str_replace(., "Human homology type","chicken vs human")
p=ggplot(df %>% filter(value == "eGene"), aes(x=Type, y=phastcons, fill=value)) + 
    #geom_violin(trim=T, fill="gray")+
    geom_boxplot(width=0.5)+
    #facet_wrap(~type, ncol=1)+
    scale_fill_manual(values=c("#F4A582", "#92C5DE"))+
    xlab("")+
    ylab("phastCons")+
    theme_classic(base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=8, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("egene_noneGene_phastCons.with_mammals.eGenes.pdf", p, width=2, height=2)



tissue_df = fread("Tissues_for_eqtl_mapping.txt") %>% as.data.frame
names(tissue_df)=c("tissue", "sample_size", "colors")
tissues=tissue_df$tissue %>% as.vector
gene_df = data.frame()
for (t in tissues){
	tmpdf = fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl_aFC.txt")) %>% as.data.frame %>% mutate(Tissue = t)
    gene_df = rbind(gene_df, tmpdf)
}
merge(df, gene_df, by.x="Gene", by.y="phenotype_id") %>% select(Gene, Type, value, log2_aFC) %>% mutate(log2_aFC=abs(log2_aFC)) -> plot_df

p=ggplot(plot_df %>% filter(value=="eGene"), aes(x=Type, y=log2_aFC, fill=value)) + 
    #geom_violin(trim=T, fill="gray")+
    geom_boxplot(width=0.5)+
    #facet_wrap(~type, ncol=1)+
    scale_fill_manual(values=c("#F4A582", "#92C5DE"))+
    xlab("")+
    ylab(expression("|log"[2]*"aFC|"))+
    theme_classic(base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=8, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("egene_noneGene_log2afc.with_mammals.eGene.pdf", p, width=2, height=2)


tissue_df = fread("Tissues_for_eqtl_mapping.txt") %>% as.data.frame
names(tissue_df)=c("tissue", "sample_size", "colors")
tissues=tissue_df$tissue %>% as.vector
gene_df = data.frame()
for (t in tissues){
	tmpdf = fread(paste0(t,"/ChickenGTEx.", t,".cis_qtl_aFC.txt")) %>% as.data.frame %>% mutate(Tissue = t)
    gene_df = rbind(gene_df, tmpdf)
}
merge(df, gene_df, by.x="Gene", by.y="phenotype_id") %>% select(Gene, Type, value, tss_distance) %>% mutate(tss_distance=abs(tss_distance)) -> plot_df

p=ggplot(plot_df %>% filter(value=="eGene"), aes(x=Type, y=tss_distance/1000, fill=value)) + 
    #geom_violin(trim=T, fill="gray")+
    geom_boxplot(width=0.5)+
    #facet_wrap(~type, ncol=1)+
    scale_fill_manual(values=c("#F4A582", "#92C5DE"))+
    xlab("")+
    ylab(expression("Distance to TSS (kb)"))+
    theme_classic(base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=8, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("egene_noneGene_Tss.with_mammals.eGene.pdf", p, width=2, height=2)


tissue_df = fread("Tissues_for_eqtl_mapping.txt") %>% as.data.frame
names(tissue_df)=c("tissue", "sample_size", "colors")
tissues=tissue_df$tissue %>% as.vector
gene_df = data.frame()
for (t in tissues){
	tmpdf = fread(paste0(t,"/ChickenGTEx.", t,".cis_independent_qtl.txt.gz")) %>% as.data.frame %>% mutate(Tissue = t)
    gene_df = rbind(gene_df, tmpdf)
}
gene_df %>% group_by(Tissue, phenotype_id) %>% summarize(indep=max(rank)) ->gene_df
merge(df, gene_df, by.x="Gene", by.y="phenotype_id") %>% select(Gene, Type, value, indep) %>% mutate(indep=abs(indep)) -> plot_df

p=ggplot(plot_df %>% filter(value=="eGene"), aes(x=Type, y=indep, fill=value)) + 
    #geom_violin(trim=T, fill="gray")+
    geom_boxplot(width=0.5, outlier.shape=NA)+
    #facet_wrap(~type, ncol=1)+
    scale_fill_manual(values=c("#F4A582", "#92C5DE"))+
    xlab("")+
    ylab("# indep SNPs")+
    ylim(0, 5)+
    theme_classic(base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=8, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("egene_noneGene_indep_snps.with_mammals.eGene.pdf", p, width=2, height=2)



df_x = df
load("cis_h2.eQTLs.Rdata")
df$pval = as.numeric(as.character(df$pval))
df$h2 = as.numeric(as.character(df$h2))
df = subset(df, h2 > 0 & h2 < 1)
merge(df, df_x, by="Gene") -> plot_df

p=ggplot(plot_df %>% filter(value=="eGene"), aes(x=Type, y=h2, fill=value)) + 
    #geom_violin(trim=T, fill="gray")+
    geom_boxplot(width=0.5)+
    #facet_wrap(~type, ncol=1)+
    scale_fill_manual(values=c("#F4A582", "#92C5DE"))+
    xlab("")+
    ylab(expression(italic(cis)*"-h"^2))+
    theme_classic(base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10),
          axis.text.x = element_text(color="black", size=8, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("egene_noneGene_h2.with_mammals.eGene.pdf", p, width=2, height=2)
