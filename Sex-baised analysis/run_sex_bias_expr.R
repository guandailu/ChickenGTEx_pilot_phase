#bash: module load R/4.0.2
library(tidyverse)
library(data.table)
library("DESeq2")
library("limma")
counts = fread("../../26_eQTL_mapping/03_phenotypes/ChickenGTEx.RawCounts.gct")
counts=as.data.frame(counts)
names(counts)[1]="gene_id"
sex_samples=as.data.frame(fread("sex_deg.samples.list", header=T, sep="\t"))
sex_samples %>% filter(Sex == "Male" | Sex == "Female") %>% group_by(Tissue_s,Sex) %>% tally() %>% filter(n>=10) %>% as.data.frame() -> tissues
subset(tissues, Tissue_s %in% names(subset(table(tissues$Tissue_s),table(tissues$Tissue_s)==2))) -> tissues

for (t in unique(tissues$Tissue_s)){
cat("Analyzing ", t, "\t")
samples=subset(sex_samples, Tissue_s==t & (Sex == "Female" | Sex == "Male"))[,"BioSample"]
counts %>% select(.,c("gene_id", samples)) -> tissue_counts
rownames(tissue_counts) = tissue_counts$gene_id
tissue_counts$gene_id=NULL
tissue_samples = filter(sex_samples, BioSample %in% samples)
tissue_samples=select(tissue_samples, c("BioSample", "Sex", "BioProject", "Release_Year", "Age/days", "Breed", "Instrument", "LibraryLayout","LibrarySelection1"))
tissue_samples=as.data.frame(tissue_samples)
rownames(tissue_samples)=tissue_samples$BioSample
tissue_samples$BioSample=NULL
tissue_samples$Sex <- factor(tissue_samples$Sex)
tissue_samples$BioProject <- factor(tissue_samples$BioProject)

all(rownames(tissue_samples) %in% colnames(tissue_counts))
names(tissue_samples)=c("Sex", "BioProject", "Year", "Age", "Breed", "Instrument", "LibraryLayout", "LibrarySelection")
tissue_samples$covar=paste(tissue_samples$BioProject,tissue_samples$Age,tissue_samples$Breed, sep="_")
tissue_samples=select(tissue_samples, Sex, covar)
tissue_samples$covar=str_replace_all(tissue_samples$covar, " ",".")
tissue_samples$covar <- as.factor(tissue_samples$covar)
dds <- DESeqDataSetFromMatrix(countData = tissue_counts,
                              colData = tissue_samples,
                              design = ~ covar + Sex)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Sex", "Male", "Female"), alpha =0.05)
res=na.omit(res)
resOrdered <- as.data.frame(res[order(res$pvalue),])
sigRes = subset(resOrdered, padj<= 0.05)
write.table(resOrdered,paste0("sex_dif_expr_analysis/",t,"_sex_DEGs.txt"), sep = "\t", quote=F, row.names=T)
write.table(sigRes,paste0("sex_dif_expr_analysis",t,"_sex_DEGs.sig.txt"), sep = "\t", quote=F, row.names=T)
}