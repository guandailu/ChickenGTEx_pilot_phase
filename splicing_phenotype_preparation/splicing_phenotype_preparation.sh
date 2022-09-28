#!/bin/bash

###########################
#
# This script is used for preparing splicing phenotypes
# Date: Sep 15, 2022
# Dailu Guan
#
###########################

set -e
usage(){
  echo "Usage: bash $0
                [ -s|--sample_tab a file with tow columns: SmapleID and Tissue (header is require)]
                [ -e|--exon_bed a bed file with exon coordinates: chr   start   end     strand  gene_id gene_name (header is require)]
                [ -g|--gene_bed a bed file with gene coordinates: chr   start   end     strand  gene_id gene_name (no header)]
                [ -c|--chr_list a text file listing chromosomes to be analyzed]
                [ -t|--leafcutter_dir a file listing tissues to be analyzed ]
                  " 1>&2
}
exit_abnormal(){
  usage
  exit 1
}

#### Parse parameters
SHORT=s:,t:,h
LONG=sample_tab:,leafcutter_dir:,help
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")
eval set -- "$OPTS"
if [[ $# -ne 0 ]]
then
    echo "Error: Inputs files not specify !!!"
    exit_abnormal
    exit 1
fi

while :
do
  case "$1" in
    -s | --sample_tab )
      sample_tab="$2"
      shift 2
      ;;
    -e | --exon_bed )
      exon_bed="$2"
      shift 2
      ;;
    -g | --gene_bed )
      gene_bed="$2"
      shift 2
      ;;
    -c | --chr_list )
      chr_list="$2"
      shift 2
      ;;
    -s | --leafcutter_dir )
      leafcutter_dir="$2"
      shift 2
      ;;
    -h | --help)
      exit_abnormal
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done


samples=( $( cat ${sample_tab} | awk '{if (NR >= 2) print $1}' | sort | uniq | tr '\n' ' ' ) )
### convert bam to junc
mkdir -p Temp
juncfile="$(mktemp Temp/juncfiles.XXXXXXXXXX)"
for sample in ${samples[@]} ### list all smaples here
do
    echo Converting $sample to $sample.junc
        ${leafcutter_dir}/scripts/bam2junc.sh ${sample}.bam Temp/${sample}.junc  # bam files must match the pattern with sample name + ".bam"
        echo Temp/${sample}.junc >> ${juncfile}
done


### cluster junc for all samples
python3 cluster_prepare_fastqtl.1.py \
        ${juncfile} \
        ${exon_bed} \
        ${gene_bed} \
        All \
        --min_clu_reads 30 --min_clu_ratio 0.001 --max_intron_len 500000 \
        --leafcutter_dir ${leafcutter_dir} \
        --output_dir All \
        --chr_list ${chr_list}


### Run filtering at tissue level
tissues=( $( cat ${sample_tab} | awk '{if (NR >= 2) print $2}' | sort | uniq | tr '\n' ' ' ) )
for tissue in ${tissues[@]} # list all tissues you want to analyze
do
  mkdir -p ${tissue}
  tissue_samples_list_file="$(mktemp ${tissue}/${tissue}.samples.XXXXXXXXXX)"
  cat ${sample_tab} | awk -v t=${tissue} '{if ($2 == t) print $1}' > ${tissue_samples_list_file}
  python3 cluster_prepare_fastqtl.2.py \
          ${tissue_samples_list_file} \  # a file with single coloumn listing sample name
          ${exon_bed} \
          ${gene_bed} \
          ${tissue} \
          --tissue ${tissue} \
          --num_pcs 10 \
          --leafcutter_dir ${leafcutter_dir} \
          --output_dir ${tissue}/
   sleep 30
   rm -rf ${tissue_samples_list_file}
done

### remove temporary folder
rm -rf Temp

### output splicing phenotype file: ${tissue}/${tissue}.leafcutter.sorted.bed.gz, and index file: ${tissue}/${tissue}.leafcutter.sorted.bed.gz.tbi

echo -e "Done!\n"
