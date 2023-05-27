#!/bin/bash


sample=$1
cellranger count --id=${sample}_res \
                   --transcriptome=/group/hjzhou98grp/dguan/80_single_cell/00_ref/Gallus_gallus.genome \
                   --fastqs=/group/hjzhou98grp/dguan/80_single_cell/01_raw_data/${sample} \
                   --sample=${sample} \
 		   --localcores=12 \
		   --chemistry threeprime



#cellranger count --id=${sample}_res --transcriptome=/group/hjzhou98grp/dguan/80_single_cell/00_ref/Gallus_gallus.genome --fastqs=/group/hjzhou98grp/dguan/80_single_cell/01_raw_data/${sample}/${sample} --sample=bamtofastq --localcores=12
