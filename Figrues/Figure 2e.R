library(data.table)
library(tidyverse)
### lead variants
tissues=fread("Tissues_for_eqtl_mapping.txt", header=F) %>% as.data.frame %>% pull(V1)
mydf = data.frame()
#tissues=c("Muscle")
for (e in c("eQTL", "sQTL", "lncQTL", "exQTL", "3aQTL")){
  if (e == "eQTL"){
      path="02_eQTLs"
  }else if (e == "sQTL"){
    path="03_sQTLs" 
  }else if (e == "lncQTL"){
    path="04_lncQTLs"
  }else if (e == "exQTL"){
    path="05_exQTLs"
  }else if (e == "3aQTL"){
    path="06_3aQTLs"
  }
  for (t in tissues){
    file=paste0(t, "/ChickenGTEx.", t,".cis_qtl.fdr.txt")
    if (file.exists(file)){
      df = fread(file, header=T) %>% as.data.frame %>% mutate(Tissue=t, QTL=e)
      if (e == "eQTL"){
        df = df %>% select(phenotype_id, variant_id, Tissue, QTL)
      }else if (e == "sQTL"){
        df = df %>% separate(phenotype_id, into=c("chr", "start", "end", "id", "phenotype_id"), sep=":") %>% select(phenotype_id, variant_id, Tissue, QTL)
      }else if (e == "lncQTL"){
        df = df %>% select(phenotype_id, variant_id, Tissue, QTL)
      }else if (e == "exQTL"){
        df = df %>% separate(phenotype_id, into=c("phenotype_id", "chr", "start", "end"), sep="_") %>% select(phenotype_id, variant_id, Tissue, QTL)
      }else if (e == "3aQTL"){
        df = df %>% separate(phenotype_id, into=c("trxid", "phenotype_id", "chr", "strand"), sep="\\|") %>% select(phenotype_id, variant_id, Tissue, QTL)
      }
      mydf = rbind(mydf, df) 
    }
  }
}
bim=fread("ChickenGTEx.bim", header=F) %>% as.data.frame
names(bim)=c("varchr", "variant_id", "pos0", "varposition", "alt", "ref")
merge(mydf, bim, by = "variant_id") -> mydf

gene_df = fread("genes.bed") %>% as.data.frame
names(gene_df)=c("Chr", "Start", "End", "Strand", "GeneID", "gene_name", "biotype")
gene_df = gene_df %>% select(GeneID, Chr, Start, End, Strand)
gene_df2=fread("lncRNA.bed", header=T) %>% as.data.frame
gene_df = rbind(gene_df, gene_df2) %>% distinct
#fwrite(gene_df, "mol_phenotype_id.txtx", row.names=F, sep="\t", quote=F)
merge(mydf, gene_df, by.x="phenotype_id", by.y="GeneID") -> mydf

subset(mydf, Strand == "+") %>% mutate(type = case_when(varposition < Start ~ "upstream",
                                                        varposition > End ~ "downstream",
                                                        varposition > Start & varposition < End ~ "body"
                                                        )) -> df1
subset(mydf, Strand == "-") %>% mutate(type = case_when(varposition > End ~ "upstream",
                                                        varposition >= Start & varposition <= End ~ "body",
                                                        varposition < Start ~ "downstream")) -> df2
rbind(df1, df2) -> mydf2
mydf2 %>% filter(type == "upstream") -> df1
subset(df1, Strand == "+") %>% mutate(win= ceiling((varposition - (Start - 1000000)) / (1000000/30)) ) -> df1.1
subset(df1, Strand == "-") %>% mutate(win= ceiling(((End + 1000000)-varposition) / (1000000/30)) ) -> df1.2
rbind(df1.1, df1.2) -> df1
df1 %>% distinct %>% group_by(win, QTL) %>% summarize(num=length(unique(variant_id))) -> df1

mydf2 %>% filter(type == "body") -> df2
subset(df2, Strand == "+") %>%  mutate(win=round((varposition - Start) / ((End - Start + 1) / 10))) -> df2.1
subset(df2, Strand == "-") %>%  mutate(win=round((varposition - Start ) / ((End - Start + 1) / 10))) -> df2.2
rbind(df2.1, df2.2) -> df2
df2$win = df2$win + 30
df2 %>% distinct %>% group_by(win, QTL) %>% summarize(num=length(unique(variant_id))) -> df2

mydf2 %>% filter(type == "downstream") -> df3
subset(df3, Strand == "+") %>% mutate(win= round((varposition - End) / (1000000/30) )) -> df3.1
subset(df3, Strand == "-") %>% mutate(win= round((Start - varposition ) / (1000000/30)) ) -> df3.2
rbind(df3.1, df3.2) -> df3
df3$win = df3$win + 40
df3 %>% distinct %>% group_by(win, QTL) %>% summarize(num=length(unique(variant_id))) -> df3

rbind(as.data.frame(df1), as.data.frame(df2), as.data.frame(df3)) -> plot_df
plot_df %>% group_by(QTL) %>% summarize(total=sum(num)) -> tdf
merge(plot_df, tdf, by="QTL") %>% mutate(pqtl=num / total * 100) -> plot_df
p=ggplot(plot_df) + 
  geom_bar(aes(x=as.integer(win), y=pqtl), stat="identity")+
  facet_wrap(~ factor(QTL, levels=c("eQTL", "sQTL", "lncQTL", "exQTL", "3aQTL")), nrow =1)+
  scale_x_continuous(breaks=c(30, 40), label=c("TSS", "TES"))+
  #scale_color_manual(values = c("#0091ff","#f0650e"))+
  xlab("Relative distance to associated genes")+
  ylab("Proportion of lead variants")+
  theme_classic(base_line_size = 0.4)+
  theme(legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=8, angle=90, hjust=1, vjust=1),
          axis.title.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=8),
          axis.title.y = element_text(color="black", size=10),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          strip.background = element_rect(color="black", fill="#C0C0C0", size=0.2, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("lead_snps_TSS_TES_distribution.pdf", p, width=7.5, height=2)

