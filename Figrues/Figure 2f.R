library(data.table)
library(tidyverse)

load("07_cis_h2/01_eQTLs/cis_h2.eQTLs.Rdata")
df$pval = as.numeric(as.character(df$pval))
df$h2 = as.numeric(as.character(df$h2))

tpm=fread("ChickenGTEx.TPM.txt") %>% as.data.frame 
rownames(tpm)=tpm$V1
tpm$V1=NULL

tissues=fread("Tissues_for_eqtl_mapping.txt", header=F) %>% as.data.frame %>% pull(V1)
mydf = data.frame()
for (t in tissues){
  df = fread(paste0(t,"/ChickenGTEx.", t,".cis_independent_qtl.aFC.txt"), header=T) %>% as.data.frame
  num_indep_df= df %>% group_by(phenotype_id) %>% summarize(num_indep=length(unique(variant_id))) %>% mutate(Tissue = t)
  merge(num_indep_df, h2_df, by.x=c("phenotype_id", "Tissue"), by.y=c("Gene", "Tissue")) -> num_indep_df
  sample_df = fread(paste0(t,"/", t,".samples.list"), header=F) %>% as.data.frame
  names(sample_df)="sampleID"
  tpm %>% select(all_of(sample_df$sampleID)) -> expr_df
  apply(expr_df, 1, median) %>% as.data.frame -> expr_df
  names(expr_df)="TPM"
  expr_df$phenotype_id = rownames(expr_df)
  merge(num_indep_df, expr_df, by="phenotype_id") -> num_indep_df
  mydf = rbind(mydf, num_indep_df)
}
mydf$num_indep=ifelse(mydf$num_indep <= 3, mydf$num_indep, ">=4")
mydf$log2TPM = log2(mydf$TPM + 0.1)
mydf = na.omit(mydf)
mydf %>% select(num_indep, TPM, log2TPM, h2) -> mydf

nonegenes=fread("non-eGenes.list", header=F) %>% as.data.frame %>% pull(V1)
tpm[nonegenes,] -> nonegenes_df
apply(nonegenes_df, 1, median) %>% as.data.frame -> nonegenes_df
names(nonegenes_df)="TPM"
nonegenes_df$Gene = rownames(nonegenes_df)
nonegenes_df$log2TPM = log2(nonegenes_df$TPM + 0.1)
subset(h2_df, Gene %in% nonegenes) -> nonegenes_h2_df
merge(nonegenes_h2_df, nonegenes_df, by="Gene") %>% mutate(num_indep=0) %>% select(num_indep, TPM, log2TPM, h2) -> neg_df

rbind(mydf, neg_df) -> mydf
ylim.prim <- c(0, 1)
ylim.sec <- c(ceiling(min(mydf$log2TPM))-1, ceiling(max(mydf$log2TPM)))
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
mydf$num_indep = factor(mydf$num_indep, levels=c("1", "2", "3", ">=4", "0"))
mydf$plot_tpm = a + mydf$log2TPM * b 

mydf %>% gather(Type, value, h2, plot_tpm) -> plot_df
p=ggplot()+
    geom_violin(plot_df, mapping=aes(x=num_indep , y=value, fill=Type), width=0.8, position=position_dodge(width=0.9))+
    geom_boxplot(plot_df, mapping=aes(x=num_indep , y=value, fill=Type), alpha=0.8, width=0.1, outlier.shape=NA, position=position_dodge(width=0.9))+
    #facet_wrap(~ann_type)+
    xlab("# independent cis-eQTLs")+
    #ylab(expression(italic(cis-h^2)))+
    scale_fill_manual(values=c("#4393C3", "#D6604D"))+
    scale_y_continuous(name = expression(italic(cis-h^2)), limits=c(0, 1), sec.axis = sec_axis(trans=~(.-a)/b , name=expression("log"[2]*"(Median TPM+0.1)")))+
    ggplot2::theme_classic(base_size = 11, base_line_size = 0.5)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.text.x = element_text(color="black", size=11),
          axis.text.y.left = element_text(color="black", size=11),
          axis.ticks.y.right=element_line(color="red"),
          axis.line.y.right=element_line(color="red"),
          axis.title.y.right=element_text(color="red", size=11),
          axis.text.y.right = element_text(color="red", size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("Num_indep_h2_expr.with_nonegene.pdf", p, width=4.5, height=3)

