library(tidyverse)
ARGS <- commandArgs(trailingOnly = TRUE)
file_perm = ARGS[1] # permutation result file
file_output = ARGS[2] # output file name
fdr_thresholds = as.numeric(ARGS[3])
# fdr_thresholds = 0.05

#Read data
eqtl = read.table(gzfile(file_perm), sep = "\t", header= T)
#eqtl = eqtl[which(!is.na(eqtl$pval_beta)),]
cat("  * Number of molecular phenotypes =", nrow(eqtl), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(eqtl$pval_beta, eqtl$pval_perm), 4), "\n")

#Run qvalue on pvalues for best signals
#Q = qvalue(eqtl$pval_beta)
#cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Alternative: Run p.adjust (BH) on pvalues for best signals
qval = p.adjust(eqtl$pval_beta, method = "BH")

#Determine significance threshold
set0 = eqtl[which(qval <= fdr_thresholds),]
set1 = eqtl[which(qval > fdr_thresholds),]
pthreshold = (sort(set1$pval_perm)[1] - sort(-1.0 * set0$pval_perm)[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")

#Calculate nominal pvalue thresholds; binominal
nthresholds = qbeta(pthreshold, eqtl$beta_shape1, eqtl$beta_shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)

eqtl$qval = qval
eqtl$pval_nominal_threshold = nthresholds
eqtl$is_eGene = (eqtl$pval_nominal < eqtl$pval_nominal_threshold) & (eqtl$qval <= fdr_thresholds)
eqtl_sig=subset(eqtl,is_eGene=="TRUE")
#Output
write.table(eqtl_sig, file_output, sep = "\t", quote=F, row.names=F)
write.table(eqtl, str_replace(file_perm,".gz",""), sep = "\t", quote=F, row.names=F)
cat("Done\n")
