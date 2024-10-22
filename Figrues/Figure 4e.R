library(tidyverse)
library(data.table)

sc_sharing = readRDS("ieqtl/pairwise_sharing.RDS")
bk_sharing = readRDS("eqtl/pairwise_sharing.RDS")

sc_r2 = NULL
for (i in 1:(nrow(sc_sharing) - 1)){
    sc_r2 = append(sc_r2, sc_sharing[i, (i+1):ncol(sc_sharing)])
}

bk_r2 = NULL
for (i in 1:(nrow(bk_sharing) - 1)){
    bk_r2 = append(bk_r2, bk_sharing[i, (i+1):ncol(bk_sharing)])
}
mydf = rbind(data.frame(Type = "eQTL", pcor = bk_r2),
             data.frame(Type = "cieQTL", pcor = sc_r2))
p=ggplot(mydf, aes(x = pcor, y=..scaled.., fill = Type)) +
    geom_density(alpha=0.5) +
    geom_vline(xintercept = median(subset(mydf, Type == "eQTL")$pcor), linetype="dashed", color = "red", linewidth=1)+
    annotate(geom="text", x= median(subset(mydf, Type == "eQTL")$pcor) + 0.1, y=1, label=round(median(subset(mydf, Type == "eQTL")$pcor), 2), color="red")+
    geom_vline(xintercept = median(subset(mydf, Type == "cieQTL")$pcor), linetype="dashed", color = "red", linewidth=1)+
    annotate(geom="text", x= median(subset(mydf, Type == "cieQTL")$pcor) + 0.1, y=1, label= round(median(subset(mydf, Type == "cieQTL")$pcor), 2), color="red")+
    xlab("MashR Spearman's rho")+
    ylab("Proportion of paired tissues/cell types")+
    scale_fill_manual(values=c("#B2182B", "#2166AC"))+
    theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.title=element_blank(),
       legend.position = c(0.8, 0.8),
       axis.text.x = element_text(color="black", size=10),
       axis.text.y = element_text(color="black", size=10),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("MashR_rho.pdf", p, width = 5, height=4)


