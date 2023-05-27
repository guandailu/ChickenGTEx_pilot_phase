library(tidyverse)
library(data.table)
library(ggsci)
library(limma)
library(pheatmap)

#### Input data
DF=fread("chicken_gtex_tpm_matrix.txt")
DF=as.data.frame(DF)
rownames(DF) = DF$Gene_ID
DF$Gene_ID = NULL
DF=as.data.frame(t(DF))

#### tissues with size >= 10
Tissue_s=read.table("samples_tissue_specificity_info.txt", header=T, sep="\t")
tissue_info=as.vector(subset(as.data.frame(table(Tissue_s$Tissue_s)),Freq>=10)[,1])
avail_tissues=as.vector(read.table("tissues_for_specificity_analysis.txt", header=F, sep ="\t")$V1)
tissues=tissue_info[tissue_info %in% avail_tissues]

#### Metadata
chicken_meta=as.data.frame(fread("sex_deg.samples.list", header=T, sep="\t"))
chicken_exp=DF

for (i in 1:length(tissues)){
  cat("Analyzing ", tissues[i], "\n")
  tissue_cata<-unique(chicken_meta$Tissue_categories[chicken_meta$Tissue_s==tissues[i]])
  tissue_cata
  tissue_cata_sample<-chicken_meta$BioSample[chicken_meta$Tissue_s==tissues[i]]
  Tissue_s_sample<-chicken_meta$BioSample[chicken_meta$Tissue_s==tissues[i]]
  remove_sample<-tissue_cata_sample[!tissue_cata_sample%in%Tissue_s_sample]
  remove_sample
  chicken_meta_data<-chicken_meta[!chicken_meta$BioSample%in%remove_sample,]
  Meta_data_chicken_sort<-chicken_meta_data[order(chicken_meta_data$Tissue_s),]
  Chicken_expression<-t(chicken_exp)
  Chicken_expression[1:10,1:10]
  Chicken_expression_sort<-Chicken_expression[, match(Meta_data_chicken_sort$BioSample, colnames(Chicken_expression))]
  Chicken_expression_sort[1:10,1:10]
  table(colnames(Chicken_expression_sort)==Meta_data_chicken_sort$BioSample)
  Chicken_expression_sort<-log2(Chicken_expression_sort+0.25)
  table(colnames(Chicken_expression_sort)==Meta_data_chicken_sort$BioSample)  
  tissue_index<-which(Meta_data_chicken_sort$Tissue_s==tissues[i])  
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_chicken_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  table(df)
  design <- model.matrix(~-1+factor(df))
  dim(design)
  colnames(design) <- c("Others","Tissue")
  table(design)
  rownames(design)=colnames(Chicken_expression_sort)
  design
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(Chicken_expression_sort, design) ## tissue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(Chicken_expression_sort))
  write.table(myresults,file=paste0("tissue_specificity_t_statistic/",tissues[i],".txt"),quote = F, sep="\t", row.names=T)
}
