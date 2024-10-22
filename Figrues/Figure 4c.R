library(tidyverse)
"%&%" = function(a, b) { paste0(a, b) }
lfsr_top_pairs = readRDS("eqtl/lfsr_m.s.RDS")
loglfsr_top_pairs = data.frame(lfsr_top_pairs < 0.05, geneid = gsub("(.)*,", "", rownames(lfsr_top_pairs)))
geneid_top_pairs = aggregate(loglfsr_top_pairs[[1]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$geneid
summary_gene_level_top_pairs = sapply(names(loglfsr_top_pairs)[1:28], function(x) aggregate(loglfsr_top_pairs[[x]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$x)
N_active_tissues_top_pairs = rowSums(summary_gene_level_top_pairs >= 1)
pDF = data.frame(group = c(rep("eGenes", length(N_active_tissues_top_pairs))), n_tissue = c(N_active_tissues_top_pairs))
pDF = pDF[pDF$n_tissue > 0,]
pDF %>% group_by(group) %>% mutate(num_tissues = cut(n_tissue, c(0,4,8,12,16,20,24,28))) %>% group_by(group, num_tissues) %>% summarize(n=n())  %>% mutate(freq=round(n/sum(n), 3)) -> edf

rm(list=ls())
lfsr_top_pairs = readRDS("sqtl/lfsr_m.s.RDS")
loglfsr_top_pairs = data.frame(lfsr_top_pairs < 0.05, geneid = gsub("(.)*,", "", rownames(lfsr_top_pairs)))
geneid_top_pairs = aggregate(loglfsr_top_pairs[[1]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$geneid
summary_gene_level_top_pairs = sapply(names(loglfsr_top_pairs)[1:27], function(x) aggregate(loglfsr_top_pairs[[x]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$x)
N_active_tissues_top_pairs = rowSums(summary_gene_level_top_pairs >= 1)
pDF = data.frame(group = c(rep("sGenes", length(N_active_tissues_top_pairs))), n_tissue = c(N_active_tissues_top_pairs))
pDF = pDF[pDF$n_tissue > 0,]
pDF %>% group_by(group) %>% mutate(num_tissues = cut(n_tissue, c(0,4,8,12,16,20,24,28))) %>% group_by(group, num_tissues) %>% summarize(n=n())  %>% mutate(freq=round(n/sum(n), 3)) -> sdf


lfsr_top_pairs = readRDS("lncqtl/lfsr_m.s.RDS")
loglfsr_top_pairs = data.frame(lfsr_top_pairs < 0.05, geneid = gsub("(.)*,", "", rownames(lfsr_top_pairs)))
geneid_top_pairs = aggregate(loglfsr_top_pairs[[1]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$geneid
summary_gene_level_top_pairs = sapply(names(loglfsr_top_pairs)[1:28], function(x) aggregate(loglfsr_top_pairs[[x]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$x)
N_active_tissues_top_pairs = rowSums(summary_gene_level_top_pairs >= 1)
pDF = data.frame(group = c(rep("lncGenes", length(N_active_tissues_top_pairs))), n_tissue = c(N_active_tissues_top_pairs))
pDF = pDF[pDF$n_tissue > 0,]
pDF %>% group_by(group) %>% mutate(num_tissues = cut(n_tissue, c(0,4,8,12,16,20,24,28))) %>% group_by(group, num_tissues) %>% summarize(n=n())  %>% mutate(freq=round(n/sum(n), 3)) -> ldf


lfsr_top_pairs = readRDS("exqtl/lfsr_m.s.RDS")
loglfsr_top_pairs = data.frame(lfsr_top_pairs < 0.05, geneid = gsub("(.)*,", "", rownames(lfsr_top_pairs)))
geneid_top_pairs = aggregate(loglfsr_top_pairs[[1]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$geneid
summary_gene_level_top_pairs = sapply(names(loglfsr_top_pairs)[1:28], function(x) aggregate(loglfsr_top_pairs[[x]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$x)
N_active_tissues_top_pairs = rowSums(summary_gene_level_top_pairs >= 1)
pDF = data.frame(group = c(rep("exGenes", length(N_active_tissues_top_pairs))), n_tissue = c(N_active_tissues_top_pairs))
pDF = pDF[pDF$n_tissue > 0,]
pDF %>% group_by(group) %>% mutate(num_tissues = cut(n_tissue, c(0,4,8,12,16,20,24,28))) %>% group_by(group, num_tissues) %>% summarize(n=n())  %>% mutate(freq=round(n/sum(n), 3)) -> xdf


lfsr_top_pairs = readRDS("3aqtl/lfsr_m.s.RDS")
loglfsr_top_pairs = data.frame(lfsr_top_pairs < 0.05, geneid = gsub("(.)*,", "", rownames(lfsr_top_pairs)))
geneid_top_pairs = aggregate(loglfsr_top_pairs[[1]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$geneid
summary_gene_level_top_pairs = sapply(names(loglfsr_top_pairs)[1:24], function(x) aggregate(loglfsr_top_pairs[[x]], by = list(geneid = loglfsr_top_pairs$geneid), sum)$x)
N_active_tissues_top_pairs = rowSums(summary_gene_level_top_pairs >= 1)
pDF = data.frame(group = c(rep("3aGenes", length(N_active_tissues_top_pairs))), n_tissue = c(N_active_tissues_top_pairs))
pDF = pDF[pDF$n_tissue > 0,]
pDF %>% group_by(group) %>% mutate(num_tissues = cut(n_tissue, c(0,4,8,12,16,20,24))) %>% group_by(group, num_tissues) %>% summarize(n=n())  %>% mutate(freq=round(n/sum(n), 3)) -> adf

rbind(edf, sdf, ldf, xdf, adf) %>% na.omit %>% as.data.frame -> mydf
mydf$group = factor(mydf$group, levels=c("eGenes", "lncGenes", "exGenes", "sGenes", "3aGenes"))
p=ggplot(mydf) +
    geom_bar(aes(x=num_tissues, y = freq, fill = group), position="dodge", stat="identity")+
    scale_fill_manual(values=c("#B2182B" , "#D6604D", "#F4A582", "#92C5DE", "#4393C3"))+
    xlab("Tissues with LFSR < 0.05")+
    ylab("Fraction of ePhenotypes")+
    theme_classic(base_size = 18, base_line_size = 1)+
    theme(legend.position = c(0.6, 0.6), legend.title=element_blank(),
         axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
         axis.text.y = element_text(color="black", size=15),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("Tissues_sharing.lfsr_0.05.all_sphenotypes.pdf", p, width = 6, height = 4)
