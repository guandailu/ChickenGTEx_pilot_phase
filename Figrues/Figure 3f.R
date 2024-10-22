library(tidyverse)
library(data.table)

mydf = fread("functional_enrichment.res.txt", header=T) %>% as.data.frame
mydf %>% group_by(Ann, QTL) %>% summarize(mFC = mean(FC), lFC = sd(FC), hFC = sd(FC)) -> plot_df
plot_df$mFC=log2(plot_df$mFC)
plot_df$lFC=log2(plot_df$lFC)
plot_df$hFC=log2(plot_df$hFC)
plot_df$Ann=factor(plot_df$Ann, levels=c(paste0("E", 1:15)))
plot_df$QTL=factor(plot_df$QTL, levels=c("eQTL", "sQTL", "lncQTL", "exQTL", "3aQTL"))

p=ggplot(plot_df) +
    geom_pointrange(aes(x=Ann, y = mFC, ymin=mFC-lFC, ymax=mFC+hFC, color = QTL), position=position_dodge(width=0.5))+
    scale_color_manual(values=c("#B2182B" , "#D6604D", "#F4A582", "#92C5DE", "#4393C3"))+
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    ylab(expression("Log"[2]*"(Fold Enrichment)"))+
    xlab("Chromatin states")+
    coord_flip()+
    theme_classic(base_size = 25, base_line_size = 1)+
    theme(
          axis.text.x = element_text(color="black", size=20),
          axis.text.y = element_text(color="black", size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("functional_enrichmentres.pdf", p, width=8, height=8)



library(data.table)
library(tidyverse)
library(ggbreak)

mydf = fread("functional_enrichment.res.txt", header=T) %>% as.data.frame
mydf$prop = mydf$num_var_in_RE / mydf$total_var
mydf %>% group_by(Ann,QTL) %>% summarize(m_prop=mean(prop), s_prop=sd(prop)) -> plot_df

plot_df$Ann=factor(plot_df$Ann, levels=c(paste0("E", 1:15)))
plot_df$QTL=factor(plot_df$QTL, levels=c("eQTL", "sQTL", "lncQTL", "exQTL", "3aQTL"))

p=ggplot(plot_df, aes(x=Ann, y = m_prop, fill = QTL)) +
    geom_errorbar(aes(ymin = m_prop, ymax=m_prop+s_prop), color="black",position = position_dodge(.9), width = 0.2)+
    geom_col(position = position_dodge())+
    scale_fill_manual(values=c("#B2182B" , "#D6604D", "#F4A582", "#92C5DE", "#4393C3"))+
    xlab("Chromatin states")+
    ylab("Proportion of variants")+
    scale_y_break(c(0.2,0.5), space=0.2, ticklabels=c(0.5, 0.6, 0.7))+
    coord_flip()+
    theme_classic(base_size = 25, base_line_size = 1)+
    theme(
          axis.text.x = element_text(color="black", size=20, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("functional_enrichment.variant_proportion.barplot.pdf", p, width=6, height=8)

