import pandas as pd
sample_file = config["sample_file"]
samples_df = pd.read_csv(sample_file).set_index("Samples", drop=False)
SAMPLES = samples_df.index.values
sample = SAMPLES
GTFFILE = config["gtf"]
REFFA = config["reffa"]
STARIDX = config["staridx"]
CHRNAME = STARIDX+"chrName.txt"
DBSNP = config["dbsnp"]
MAP = config["map"]
intervals = config["intervals"]
interval_df = pd.read_csv(intervals).set_index("Interval", drop=False)
INTERVALS = [str(i) for i in interval_df.index.values]
interval = INTERVALS
MAPPED = expand("03_aligned/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES)
MAPRE = expand("03_aligned/{sample}_Log.final.out", sample = SAMPLES)
BAI = expand("03_aligned/{sample}_Aligned.sortedByCoord.out.bam.bai",sample =SAMPLES)
STTGTF = expand("04_stringtie_quant/{sample}.gtf", sample = SAMPLES)
STTAE = expand("04_stringtie_quant/{sample}_stringtie.tsv", sample = SAMPLES)
FTCOUNTS = expand("05_featureCounts_quant/{sample}_featureCounts_quant.txt",sample = SAMPLES)
RGBAM = expand("06_addrg/{sample}_rg.bam", sample = SAMPLES)
RGBAI = expand("06_addrg/{sample}_rg.bai", sample = SAMPLES)
JUNC = expand("06_addrg/{sample}_rg.bam.junc",sample = SAMPLES)
MKBAM = expand("07_mkdup/{sample}_mkdup.bam", sample = SAMPLES)
MKBAI = expand("07_mkdup/{sample}_mkdup.bai", sample = SAMPLES)
CIGARBAM = expand("08_cigar/{sample}_cigar.bam", sample = SAMPLES)
CIGARBAI = expand("08_cigar/{sample}_cigar.bai", sample = SAMPLES)
BQSR = expand("09_bqsr/{sample}_bqsr.bam", sample = SAMPLES)
BAMIDX = expand("09_bqsr/{sample}_bqsr.bai", sample = SAMPLES)
BAMMD5 = expand("09_bqsr/{sample}_bqsr.bam.md5", sample = SAMPLES)
GVCF = expand("10_indiv_gvcf/{sample}_SNPs_called.g.vcf.gz", sample = SAMPLES)
GVCFIDX = expand("10_indiv_gvcf/{sample}_SNPs_called.g.vcf.gz.tbi", sample = SAMPLES)
iVCF = expand("11_RNAseq_GenomicsDB/Sheep_GTEx_chr{interval}", interval = INTERVALS)
VCF = expand("12_GenotypeGVCFs_vcf/Sheep_GTEx_chr{interval}.vcf.gz", interval = INTERVALS)
RVCF = "13_Filtration_vcf/Sheep_GTEx.vcf.gz"
FVCF = "13_Filtration_vcf/Sheep_GTEx_filtrated.vcf.gz"
SNPs = "14_select_vcf/Sheep_GTEx_filtrated.SNPs.vcf.gz"
SNPidx = "14_select_vcf/Sheep_GTEx_filtrated.SNPs.vcf.gz.tbi"
ASE = expand("15_ASE_count/{sample}_ASE_table.txt", sample = SAMPLES)
TVCF = "13_Filtration_vcf/Sheep_GTEx.vcf.gz.tbi"


if config["paired_reads"] is True:
	ruleorder: all > build_star_index > trim_paired_reads > paired_star_align > idxbam > stringtie_quant > featureCounts_quant > addrg > bam2junc > markduplicates > splitncigarreads > baserecalibrator > applybqsr > haplotypcaller2gvcf > build_database > genotyping > bcftools > IndexFeatureFile >  HardFilter > selectSNP > ASEcount
rule all:
    input:
    	MAPPED, MAPRE, BAI, STTGTF, STTAE, FTCOUNTS, RGBAM, RGBAI, JUNC,  MKBAM, MKBAI, CIGARBAM, CIGARBAI, BQSR, BAMIDX, BAMMD5, GVCF, GVCFIDX, iVCF, VCF, RVCF, FVCF, SNPs, SNPidx, ASE, TVCF

rule build_star_index:
    input:
        gtf = GTFFILE,
        fa = REFFA
    output:
        CHRNAME
    params:
        starref = STARIDX
    conda:
        "envs/star.yaml"
    threads: 4
    message: "Start building STAR index..."
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.starref} --runThreadN {threads} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf}"

if config["paired_reads"] is True:
    rule trim_paired_reads:
    	input:
            R1 = "PRJEB19199/{sample}_R1.fastq.gz",
            R2 = "PRJEB19199/{sample}_R2.fastq.gz"
    	output:
            TR1 = "01_trimmed/{sample}_val_1.fq.gz",
            TR2 = "01_trimmed/{sample}_val_2.fq.gz"
    	params:
            basename = "{sample}",
            outdir = "01_trimmed"
    	log:
	    "{sample}.log"
    	conda:
            "envs/trim_galore.yaml"
    	threads: 12
    	message: "Start trimming PAIRED reads..."
    	shell:
            """
            trim_galore --paired --gzip --trim-n --cores {threads} --length 30 --clip_R1 3 --clip_R2 3 --three_prime_clip_R1 3 --three_prime_clip_R2 3 -o {params.outdir} --basename {params.basename} {input.R1} {input.R2}
            """
if config["paired_reads"] is True:
    rule paired_star_align:
        input:
            R1 = "01_trimmed/{sample}_val_1.fq.gz",
            R2 = "01_trimmed/{sample}_val_2.fq.gz",
            gtf = GTFFILE,
            staridx = STARIDX,
        output:
            bam = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam",
            logout = "03_aligned/{sample}_Log.final.out"
        params:
            prefix = "03_aligned/{sample}"+"_"
        conda:
            "envs/star.yaml"
        threads: 12
        message: "Start runing alignment..."
        shell:
            """
            STAR --runThreadN {threads} --genomeDir {input.staridx} --sjdbGTFfile {input.gtf} --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.prefix} --quantMode GeneCounts --chimSegmentMin 10 --chimOutType Junctions --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --readFilesCommand zcat --outFilterMismatchNmax 3
            """

rule idxbam:
    input:
        bam = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bai = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam.bai"
    threads: 8
    conda:
        "envs/samtools.yaml"
    message: "Start indexing bam..."
    shell:
        "samtools index -@ {threads} {input}"

rule stringtie_quant:
    input:
        bam = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = GTFFILE
    output:
        samplegtf = "04_stringtie_quant/{sample}.gtf",
        ae =  "04_stringtie_quant/{sample}_stringtie.tsv"
    conda:
        "envs/stringtie.yaml"
    threads: 12
    message: "Start quantifying gene expression by StringTie..."
    shell:
        "stringtie -e -B -p {threads} -G {input.gtf} -o {output.samplegtf} -A {output.ae} {input.bam}"

rule featureCounts_quant:
    input:
        bam = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = GTFFILE
    output:
        re = "05_featureCounts_quant/{sample}_featureCounts_quant.txt",
        sm = "05_featureCounts_quant/{sample}_featureCounts_quant.txt.summary"
    conda:
        'envs/featureCounts.yaml'
    threads: 12
    shell:
        "featureCounts -f -t exon -g gene_id -a {input.gtf} -T {threads} -o {output.re} {input.bam}"

rule addrg:
    input:
        ibam = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        ibai = "03_aligned/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        obam = "06_addrg/{sample}_rg.bam",
        obai = "06_addrg/{sample}_rg.bai"
    params:
        sm = "{sample}"
    threads: 4
    conda:
        "envs/gatk4.yaml"
    message: "Start adding read groups..."
    shell:
        """
        gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/tmp" AddOrReplaceReadGroups -I {input.ibam} -O {output.obam} -LB {params.sm} -PL ILLUMINA -PU {params.sm} -SM {params.sm} -ID {params.sm} -FO {params.sm} --CREATE_INDEX true
        """

rule bam2junc:
    input:
        bam = "06_addrg/{sample}_rg.bam"
    output:
        junc = "06_addrg/{sample}_rg.bam.junc"
    threads: 4
    message: "Converting bamfile to bamfile.junc"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        /home/sheep/feng/r1-test/leafcutter/leafcutter-master/scripts/bam2junc.sh {input.bam}  {output.junc}
        """

rule markduplicates:
    input:
        ibam = "06_addrg/{sample}_rg.bam",
        ibai = "06_addrg/{sample}_rg.bai"
    output:
        obam = "07_mkdup/{sample}_mkdup.bam",
        obai = "07_mkdup/{sample}_mkdup.bai",
        metrics = "07_mkdup/{sample}_metrics.txt"
    threads: 4
    conda:
        "envs/gatk4.yaml"
    message: "Start marking duplicates..."
    shell:
        """
        gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/tmp -XX:-UseGCOverheadLimit" MarkDuplicates --spark-runner LOCAL -I {input.ibam} -O {output.obam} -M {output.metrics} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.0025
        """

rule splitncigarreads:
    input:
        ibam = "07_mkdup/{sample}_mkdup.bam",
        ibai = "07_mkdup/{sample}_mkdup.bai",
        ref = REFFA    
    output:
        obam = "08_cigar/{sample}_cigar.bam",
        obai = "08_cigar/{sample}_cigar.bai"
    params:
        chr = config["genome_chr"]
    threads: 12
    conda:
        "envs/gatk4.yaml"
    message: "Start running SplitNCigarReads..."
    shell:
        """
        gatk --java-options "-Xmx64G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=/tmp" SplitNCigarReads --spark-runner LOCAL -I {input.ibam} -R {input.ref} -O {output.obam} --create-output-bam-index true --max-reads-in-memory 1000
        """

rule baserecalibrator:
    input:
        ibam = "08_cigar/{sample}_cigar.bam",
        ibai = "08_cigar/{sample}_cigar.bai",
        ref = REFFA
    output:
        bqsr = "09_bqsr/{sample}_bqsr.table"
    params:
        dbsnp = DBSNP
    threads: 12
    conda:
        "envs/gatk4.yaml"
    message: "Start running BaseRecalibrator..."
    shell:
        """
        gatk --java-options "-Xmx64G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=/tmp" BaseRecalibrator --spark-runner LOCAL -I {input.ibam} --known-sites {params.dbsnp} -O {output} -R {input.ref} --tmp-dir /tmp
        """

rule applybqsr:
    input:
        ibam = "08_cigar/{sample}_cigar.bam",
        ibai = "08_cigar/{sample}_cigar.bai",
        bqsr = "09_bqsr/{sample}_bqsr.table",
        ref = REFFA
    output:
        obam = "09_bqsr/{sample}_bqsr.bam",
        obai = "09_bqsr/{sample}_bqsr.bai",
        md5 = "09_bqsr/{sample}_bqsr.bam.md5"
    params:
        dbsnp = DBSNP  # vcf file should be indexed
    threads: 12
    conda:
        "envs/gatk4.yaml"
    message: "Start running ApplyBQSR..."
    shell:
        """
        gatk --java-options "-Xmx64G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=/tmp" ApplyBQSR --spark-runner LOCAL -I {input.ibam} --bqsr-recal-file {input.bqsr} -O {output.obam} --create-output-bam-index true --create-output-bam-md5 true --tmp-dir /tmp
        """

rule haplotypcaller2gvcf:
    input:
        ibam ="08_cigar/{sample}_cigar.bam",
        ref = REFFA
    output:
        gvcf = "10_indiv_gvcf/{sample}_SNPs_called.g.vcf.gz",
        vcfidx = "10_indiv_gvcf/{sample}_SNPs_called.g.vcf.gz.tbi"
    threads: 12
    conda:
        "envs/gatk4.yaml"
    message: "Start generating GVCFs..."
    shell:
        """
        gatk --java-options "-Xmx64G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=/tmp" HaplotypeCaller --spark-runner LOCAL --native-pair-hmm-threads {threads} -ERC GVCF -I {input.ibam} -O {output.gvcf}  -R {input.ref} --create-output-variant-index true --all-site-pls true --dont-use-soft-clipped-bases true --output-mode EMIT_ALL_CONFIDENT_SITES --standard-min-confidence-threshold-for-calling 0
        """

rule build_database:
    input:
        gvcfs=expand("10_indiv_gvcf/{sample}_SNPs_called.g.vcf.gz", sample=SAMPLES),
    	map = MAP,
        ref = REFFA
    output:
        ourdir = directory("11_RNAseq_GenomicsDB/Sheep_GTEx_chr{interval}")
    params:
    	path = "11_RNAseq_GenomicsDB/Sheep_GTEx_chr{interval}",
        interval = "chr{interval}"
    conda:
        "envs/gatk4.yaml"
    threads: 12
    shell:
        """
        gatk --java-options "-Xmx24G -XX:ParallelGCThreads=12" GenomicsDBImport -R {input.ref} --sample-name-map {input.map} --genomicsdb-workspace-path {params.path} -L {params.interval} --reader-threads 12 --batch-size 50 --tmp-dir /tmp
        """

rule genotyping:
    input:
        ref = REFFA,
        db = "11_RNAseq_GenomicsDB/Sheep_GTEx_chr{interval}"
    output:
        vcf = "12_GenotypeGVCFs_vcf/Sheep_GTEx_chr{interval}.vcf.gz"
    conda:
        "envs/gatk4.yaml"
    threads: 12
    shell:
        """
        gatk --java-options "-Xmx24G -XX:ParallelGCThreads=12" GenotypeGVCFs -R {input.ref} -V gendb://{input.db} -O {output.vcf} --tmp-dir /tmp
        """

rule bcftools:
    input:
        vcfs =expand("12_GenotypeGVCFs_vcf/Sheep_GTEx_chr{interval}.vcf.gz",interval = INTERVALS),
        ref = REFFA
    output:
        rvcf = "13_Filtration_vcf/Sheep_GTEx.vcf.gz"
    conda:
        "envs/bcftools.yaml"
    threads: 12
    shell:
        """
        bcftools concat --threads {threads} -O z -o {output.rvcf} {input.vcfs}
        """
rule IndexFeatureFile:
    input:
        rvcf = "13_Filtration_vcf/Sheep_GTEx.vcf.gz"
    output:
        tvcf = "13_Filtration_vcf/Sheep_GTEx.vcf.gz.tbi"
    conda:
        "envs/gatk4.yaml"
    threads: 12
    shell:
        """
        gatk --java-options "-Xmx24G -XX:ParallelGCThreads=12" IndexFeatureFile -I {input.rvcf}
        """

rule HardFilter:
    input:
        rvcf = "13_Filtration_vcf/Sheep_GTEx.vcf.gz",
        ref = REFFA
    output:
        fvcf = "13_Filtration_vcf/Sheep_GTEx_filtrated.vcf.gz"
    conda:
        "envs/gatk4.yaml"
    threads: 12
    shell:
        """
        gatk --java-options "-Xmx24G -XX:ParallelGCThreads=12" VariantFiltration -R {input.ref} -V {input.rvcf} -O {output.fvcf} -window 35 -cluster 3 --filter-name one --filter-expression "FS>30.0" --filter-name two --filter-expression "QD<2.0" --tmp-dir /tmp
        """

rule selectSNP:
    input:
        vcf = "13_Filtration_vcf/Sheep_GTEx_filtrated.vcf.gz",
        ref = REFFA
    output:
        ovcf = "14_select_vcf/Sheep_GTEx_filtrated.SNPs.vcf.gz",
        idx = "14_select_vcf/Sheep_GTEx_filtrated.SNPs.vcf.gz.tbi"
    conda:
        "envs/gatk4.yaml"
    threads: 12
    shell:
        """
        gatk --java-options "-Xmx24G -XX:ParallelGCThreads=12" SelectVariants -R {input.ref} -V {input.vcf} -O {output.ovcf} --select-type-to-include SNP --exclude-filtered --create-output-variant-index true --tmp-dir /tmp
        """

rule ASEcount:
    input:
        bam = "06_addrg/{sample}_rg.bam",
        bai = "06_addrg/{sample}_rg.bai",
        variant = "14_select_vcf/Sheep_GTEx_filtrated.SNPs.vcf.gz",
        ref = REFFA
    output:
        ase = "15_ASE_count/{sample}_ASE_table.txt"
    threads: 12
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "-Xmx64G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=/tmp" ASEReadCounter -R {input.ref} -I {input.bam} -V {input.variant} -min-depth 10 -O {output.ase}
        """
