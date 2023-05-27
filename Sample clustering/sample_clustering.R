#### bash: module load R/4.1.0
setwd("/group/zhougrp2/dguan/07_stringtie_quant/gene_atlas/")
library(tidyverse)
library(data.table)
library("dendextend")
library("RColorBrewer")
library("maps")
library(ggtree)

#### Input data
files=list.files("../", "*_stringtie.tsv")
files=paste0("../", files)
data_list = lapply(files, fread)
data_list = lapply(data_list, function(x)x[,c("Gene ID", "TPM")])
data_list = lapply(data_list, setNames, nm=c("GeneID", "TPM"))
data_list %>% reduce(inner_join, by = "GeneID")  -> df
df_names = c("GeneID",str_replace_all(files, "../|_stringtie.tsv",""))
names(df)= df_names
df = as.data.frame(df)
dat=fread("../gene_tpm_all_samples.tsv")
dat=as.data.frame(dat)
DF=inner_join(dat, df, by = c("Gene_ID" = "GeneID"))
dim(DF)

#### Keep valid samples
samples=fread("samples_info.txt")
samples=as.data.frame(samples)
DF=select(DF, "Gene_ID", as.vector(samples$BioSample))
fwrite(DF, "chicken_gtex_tpm_matrix.txt", sep="\t")

                   
######################  Clustering using all genes    ######################                            
DF=fread("chicken_gtex_tpm_matrix.txt")
samples=fread("samples_info.txt")
samples=as.data.frame(samples)
#### Colors used in this study
tissue_colors <- read.table("Chicken_GTEx_tissue_color_setting.txt", header = T, sep = "\t", comment.char = "@")
colors <- as.vector(tissue_colors$color)
names(colors) <- as.vector(tissue_colors[match(colors, tissue_colors$color),1])

#### Log transformation of data frame
DF=as.data.frame(DF)
rownames(DF) = DF$Gene_ID
genes = as.vector(DF$Gene_ID)
DF$Gene_ID = NULL
DF=apply(log(DF+0.25), MARGIN = 2, FUN = scale)
rownames(DF) = genes
DF_sd = apply(DF, MARGIN=1, sd)
DF <- DF[order(DF_sd,decreasing=T),]

#### Tree clustering using pearson correlation of samples
res <- cor(DF[1:5000,],method="pearson")
dist <- as.dist(1 - res)
tree=hclust(dist, method="complete")
save(tree, file = "Chicken_GTEx_all_genes_expr_tree.Rdata")
samples$colors <- colors[match(samples$Tissue_categories,names(colors))]
ann_colors=samples[match(tree$labels, samples$BioSample),"colors"]

#load("Chicken_GTEx_all_genes_expr_tree.Rdata")
#### Tree viz
p=ggtree(tree, layout="fan", open.angle=90)+
    theme(legend.position='none') +
        geom_tippoint(color=ann_colors, shape=124, size=4)

pdf("Chicken_GTEx_all_genes_expr_tree.pdf", width=8, height=8)
print(p)
dev.off()