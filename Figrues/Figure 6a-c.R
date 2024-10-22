library(data.table)
library(tidyverse)
library(rstatix)

system("mkdir -p output")
twas_df = rbind(twas_pig, twas_cattle, twas_human)
tissues=fread("Tissues_eqtl_mapping.txt") %>% pull(V1)
res_df = data.frame()
for (ts in tissues){
for (t in unique(twas_df$Trait1)){
  subset(twas_df, Tissue1 == ts & Trait1 == t) -> df_twas
  df_twas$fdr <- p.adjust(df_twas$p, method="fdr")
  subset(df_twas, fdr < 0.1) -> subdf
  n_sig=nrow(subdf)
  if (n_sig > 0){
    res_df = rbind(res_df, subdf)
    source("manhattan_plot.R")
    df=df_twas
    df %>% group_by(Species2) %>% mutate(x_num=row_number()) %>% as.data.frame -> df
    df %>% group_by(Species2) %>% summarise(chr_len=max(x_num)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(df, ., by=c("Species2"="Species2")) %>% arrange(Species2, x_num) %>% mutate(BPcum=x_num+tot) -> df
    df_center = df %>% group_by(Species2) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    D <- transform(df, chr = as.numeric(as.factor(df$Tissue1)))
    chrlabs=as.vector(unique(D$Tissue1))
    fdr_threshold=D %>% filter(fdr > 0.1) %>% slice(which.min(p)) %>% pull(p)
    D$sign = ifelse(D$fdr < 0.1 & D$cor > 0, "positive", ifelse(D$fdr < 0.1 & D$cor < 0, "negative", "not"))
    shapes=c("+", "-")
    names(shapes)=c("positive", "negative")
    p=ggplot(D)+
      geom_segment(mapping=aes(x=BPcum, y=0, xend=BPcum, yend=-log10(p), color=Species2), size=0.4) +
      geom_point(mapping=aes(x = BPcum, y = -log10(p), color=Species2), size=1, shape = 19)+
      geom_point(D %>% subset(fdr< 0.1), mapping=aes(x =BPcum, y = -log10(p)), color="red", shape=19, size=4)+
      geom_point(D %>% subset(fdr< 0.1), mapping=aes(x = BPcum, y = -log10(p), shape=sign), size=4, color = "white")+
      #geom_hline(yintercept=-log10(9.11e-3), linetype="dashed", color = "blue4")+
      geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color = "red")+
      scale_color_manual(values=c("#404040", "#A0A0A0", "#404040", "#A0A0A0"))+
      scale_shape_manual(values=c(shapes))+
      xlab("")+
      ylab(expression("-log"[10]*italic(P)))+
      scale_x_continuous(breaks = c(df_center$center), label = c(df_center$Species2))+
      theme_classic(base_size = 15, base_line_size = 0.6)+
      theme(legend.title=element_blank(), legend.position="none",
            axis.text.x = element_text(color="black", size=10, angle=45, hjust=1),
            axis.text.y = element_text(color="black", size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(2, 1, 2, 2))
    ggsave(paste0("chicken_", ts,"_", t,"_three_species_traits.pdf"), p, width=3, height=2.5)
    pdf(paste0("chicken_", ts,"_", t,"_three_species_traits.qqplot.pdf"), width=3, height=4)
    library(qqman)
    qq(D$p)                                         
    dev.off() 
  }
}
} 


x=257
df = fread("chicken_mammal_traits_TWAS.txt") %>% as.data.frame

for (x in 1:nrow(df)){
Species1=df[x,"Species1"] %>% as.character
Tissue1=df[x,"Tissue1"] %>% as.character
Trait1=df[x,"Trait1"] %>% as.character
Species2=df[x,"Species2"] %>% as.character
Tissue2=df[x,"Tissue2"] %>% as.character
Trait2=df[x,"Trait2"] %>% as.character

twas1=fread(paste0(Species1, "TWAS/", Tissue1, "/", Tissue1, ".", Trait1, ".csv")) %>% as.data.frame
twas1 %>% select(gene, zscore) -> twas1
twas2=fread(paste0(Species2, "TWAS/", Tissue2, "/", Tissue2, ".", Trait2, ".csv")) %>% as.data.frame

if (Species2 == "human" | Species2 == "Human"){
    twas2  %>% separate(., gene, into=c("gene", "version"), sep="\\.") %>% select(gene, zscore) -> twas2
}else{
   twas2  %>%  select(gene, zscore) -> twas2
}
if (Species2 == "Cattle"){
  orthologs = fread("chicken_cattle_pair", sep="\t", header=T) %>% as.data.frame
}else if (Species2 == "Pig"){
  orthologs = fread("chicken_pig_pair", sep="\t", header=T) %>% as.data.frame
}else if (Species2 == "Human"){
  orthologs = fread("chicken_human_pair", sep="\t", header=T) %>% as.data.frame
}
names(orthologs)=c(Species1, Species2)

merge(orthologs, twas1, by.x=Species1, by.y="gene") -> twas1_df
merge(orthologs, twas2, by.x=Species2, by.y="gene") -> twas2_df
merge(twas1_df, twas2_df, by=Species1) -> test_df
#test_df$log_zscore.x=abs(log2(test_df$zscore.x))
#test_df$log_zscore.y=abs(log2(test_df$zscore.y))
na.omit(test_df) -> test_df
p=ggplot(test_df, aes(x = abs(zscore.x), y = abs(zscore.y)))+
                geom_point(color="gray64", fill="lightblue3", shape=21, size=3)+
                geom_smooth(method = lm, color="dodgerblue4", size=1)+
                xlab(paste0(Species1, " |Effect size|"))+
                ylab(paste0(Species2, " |Effect size|"))+
                theme_classic(base_size = 15, base_line_size = 0.9)+
                theme(axis.text.x = element_text(color="black", size=10),
                  axis.text.y=element_text(color="black", size=10),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"))
   ggsave(paste0(Species1,"_", Tissue1, "_", Trait1, ".",  Species2,"_", Tissue2, "_", Trait2, ".", "zscore.pdf"), p, width=4, height=3.5)
}

