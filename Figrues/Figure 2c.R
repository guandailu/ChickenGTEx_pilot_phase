library(data.table)
library(tidyverse)
df1=fread("LDBlocks.DR2.bed", header=F) %>% as.data.frame()
df2=fread("LDBlocks.NumGenes.bed", header=F) %>% as.data.frame()
df3 = fread("LDBlocks.SNPnum.bed", header=F) %>% as.data.frame()
names(df1)=c("chr", "start", "end", "DR2")
names(df2)=c("chr", "start", "end", "GeneNum")
names(df3)=c("chr", "start", "end", "SNPNum")
df=merge(merge(df1, df2, by=c("chr", "start", "end")), df3, by=c("chr", "start", "end"))
#df3$SNPNum_formed = ifelse(df3$SNPNum >= 50, 50, df3$SNPNum) 
p=ggplot(df3)+
    geom_histogram(aes(x=log2(SNPNum), y = ..count../1000), binwidth = 0.5) +
    geom_vline(xintercept = log2(median(df3$SNPNum)), linetype="dashed", color = "red", size=1.5)+
    #xlim(0, 50)+
    xlab(expression("log"[2]*"(# SNPs per LD block)"))+
    ylab("#LD blocks (×1000)")+
    ggplot2::theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=15),
          axis.text.y =element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("numSNPs_per_LD_block.pdf", p, width=4, height=3.5)    