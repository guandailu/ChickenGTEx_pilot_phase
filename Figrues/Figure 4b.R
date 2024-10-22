library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library("gplots")


tissue_colors <- read.table("Chicken_GTEx_colors.txt", header = T, sep = "\t", comment.char = "@")
colors <- as.vector(tissue_colors$color_hex)
names(colors) <- as.vector(tissue_colors[match(colors, tissue_colors$color_hex),"tissue_id"])

edf=readRDS("eqtl/pairwise_sharing.RDS")
sdf=readRDS("sqtl/pairwise_sharing.RDS")
ldf=readRDS("lncqtl/pairwise_sharing.RDS")
xdf=readRDS("exqtl/pairwise_sharing.RDS")
pdf=readRDS("3aqtl/pairwise_sharing.RDS")

rank_inx = function(df1, df2, qtl1, qtl2){
    df2_tissues = colnames(df2)
    if (qtl2 == "3aQTL"){
        subset(df1, rownames(df1) %in% df2_tissues) %>% as.data.frame %>% select(all_of(df2_tissues)) %>% as.matrix -> df1
    }
    res <- data.frame(qtl1=qtl1, qtl2=qtl2, rand.index(df1, df2))
    names(res)[3]="rank_id"
    return(res)
}
mydf = rbind(rank_inx(edf, sdf, "eQTL", "sQTL"),
    rank_inx(edf, ldf, "eQTL", "lncQTL"),
    rank_inx(edf, xdf, "eQTL", "exQTL"),
    rank_inx(edf, pdf, "eQTL", "3aQTL"),
    rank_inx(sdf, ldf, "sQTL", "lncQTL"),
    rank_inx(sdf, xdf, "sQTL", "exQTL"),
    rank_inx(sdf, pdf, "sQTL", "3aQTL"),
    rank_inx(ldf, xdf, "lncQTL", "exQTL"),
    rank_inx(ldf, pdf, "lncQTL", "3aQTL"),
    rank_inx(xdf, pdf, "exQTL", "3aQTL"))
p=ggplot(mydf, aes(qtl1, qtl2)) +
    geom_tile(aes(fill = rank_id), color='white') +
    xlab("")+
    ylab("")+
    scale_fill_gradient2(low = 'white', mid="darkorange2", high = 'darkred', space = 'Lab') +
    theme_classic(base_size = 24, base_line_size = 1)+
    theme(legend.title=element_blank(),
       axis.text.x = element_text(color="black", size=20),
       axis.text.y = element_text(color="black", size=20),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())
ggsave("multi_molQTLs_shairng.rank_idx.pdf",p, width=7,height=6)
