library(data.table)
library(tidyverse)
tissues=fread("Tissues_for_eqtl_mapping.txt") %>% pull(V1)
mydf = data.frame()

files=list.files("results_lm/", "_coloc.txt$")
for (f in 1:length(files)){
  file=paste0("results_lm/", files[f])
  if (file.info(file)[["size"]] != 0 ){
    df = fread(file, header=T) %>% as.data.frame
    df$Tissue =unlist(strsplit(files[f], split="\\."))[1]
    df$coloc_type=str_replace(unlist(strsplit(files[f], split="\\."))[2], "_coloc","")
    mydf = rbind(mydf, df)
  }
}
mydf %>% gather(., Type, value, PP.H4.abf, lead_ld) -> plot_df 

plot_df$coloc_type=factor(plot_df$coloc_type, levels=c("eQTL_sQTL", "eQTL_exQTL", "eQTL_apaQTL", "sQTL_exQTL", "sQTL_lncQTL", "sQTL_apaQTL", "lncQTL_exQTL", "lncQTL_apaQTL", "exQTL_apaQTL"))

lab_df = data.frame(xlab=c("eQTL vs sQTL", "eQTL vs exQTL", "eQTL vs 3aQTL", "sQTL vs exQTL", "sQTL vs lncQTL", "sQTL vs 3aQTL", "lncQTL vs exQTL", "lncQTL vs 3aQTL", "exQTL vs 3aQTL"), 
                coloc_type=c("eQTL_sQTL", "eQTL_exQTL", "eQTL_apaQTL", "sQTL_exQTL", "sQTL_lncQTL", "sQTL_apaQTL", "lncQTL_exQTL", "lncQTL_apaQTL", "exQTL_apaQTL"))
merge(plot_df, lab_df, by="coloc_type") -> plot_df


p=ggplot(plot_df, aes(x=value, fill=Type)) + 
    #geom_violin(trim=T, width=3, position = position_dodge(width = 1.5))+
    geom_histogram(alpha=0.4, position="identity", binwidth=0.05)+
    facet_wrap(~factor(xlab, c("eQTL vs sQTL", "eQTL vs exQTL", "eQTL vs 3aQTL", "sQTL vs exQTL", "sQTL vs lncQTL", "sQTL vs 3aQTL", "lncQTL vs exQTL", "lncQTL vs 3aQTL", "exQTL vs 3aQTL")), scales="free")+
    xlab("")+
    ylab("Counts")+
    scale_fill_manual(values=c("#B2182B", "#2166AC"))+
    scale_color_manual(values=c("#B2182B", "#2166AC"))+
    theme_classic(base_size = 10, base_line_size = 0.5)+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=8, angle=90, vjust=0.5, hjust=1),
          axis.text.y = element_text(color="black", size=8),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("mulQTL_coloc.pdf", p, width=6, height=4)
