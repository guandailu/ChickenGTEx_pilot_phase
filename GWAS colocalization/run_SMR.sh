#!/bin/bash

echo -e "
# =====================================================================
# Run SMR and HEIDI tests
#   * We accept genotype files in compressed vcf (.vcf.gz) and plink format
#   * We accept GWAS results from GCTA (including MLMA and fastGWA) and GEMMA
#   * We accept eqtl results from MatrixEQTL, tensorQTL and fastqtl
#   * eqtl results could be compressed (.gz) or not
#   * Note: eqtl_formats only accept arguments: "MatrixEQTL", "tensorQTL", "fastqtl"
# =====================================================================
"
set -e
usage(){
  echo "Usage: $0 [ -g|--genotype genotype file ]
                  [ -w|--gwas gwas results ]
                  [ -e|--eqtl_file eqtl results ]
                  [ -f|--eqtl_format eqtl software ]
                  [ -m|--eqtl_map eqtl plink prefix ]
                  [ -b|--eqtl_bed eqtl phenotypes ]
                  [ -t|--threads number of threads ]
                  [ -o|--out prefix of output file ] " 1>&2
}
exit_abnormal(){
  usage
  exit 1
}

#### Check if executable is available
smr="/home/dguan/bin/smr_v1.3.1_linux_x86_64_static"
if [[ ! -f ${smr} ]]; then echo 'ERROR: smr executable not available'; exit 1; fi

#### Parse parameters
SHORT=g:,w:,e:,f:,m:,b:,t:,o:,h
LONG=genotype:,gwas:,eqtl_file:,eqtl_format:,eqtl_map:,eqtl_bed:,threads:,out:,help
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")
eval set -- "$OPTS"
if [[ $? -ne 0 ]]
then
    echo "Error: Inputs files not specify !!!"
    exit_abnormal
    exit 1
fi
while :
do
  case "$1" in
    -g | --genotype )
      genotype="$2"
      shift 2
      ;;
    -w | --gwas )
      gwas="$2"
      shift 2
      ;;
    -e | --eqtl_file )
      eqtl_file="$2"
      shift 2
      ;;
    -f | --eqtl_format )
      eqtl_format="$2"
      shift 2
      ;;
    -m | --eqtl_map )
      eqtl_map="$2"
      shift 2
      ;;
    -b | --eqtl_bed )
      eqtl_bed="$2"
      shift 2
      ;;
    -t | --threads )
      threads="$2"
      shift 2
      ;;
    -o | --out )
      out="$2"
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


#### Test input files
echo -e "\n******** Check input files ********\n"
InputFiles=("${gwas}" "${eqtl_file}" "${eqtl_map}.bim" "${eqtl_bed}")
for file in ${InputFiles[@]};
do
    if [[ ! -f ${file} ]]
    then
       echo "Error: ${file} files not exist!!!"
       exit_abnormal
       exit 1
     fi
done

#### Create temporary directory
mkdir -p Temp

#### Preapre genotypes in plink format
echo -e "\n******** Prepare genotypes ********\n"
if [[ ${genotype} == *vcf.gz ]]
then
    module load plink/1.90

    bfile="$(mktemp Temp/genotypes.XXXXXXXXXX)"
    plink --vcf ${genotype} --chr-set 33 --keep-allele-order --threads ${threads} --double-id --make-bed --out ${bfile}
else
    bfile=${genotype}
fi

#### Prepare GWAS file
echo -e "\n******** Prepare GWAS results ********\n"

gwas_headers="$(mktemp Temp/gwas_headers.XXXXXXXXXX)"
gwas_prefix="$(mktemp Temp/gwas.XXXXXXXXXX)"
cat ${gwas} | head -n 1 | tr '\t' '\n' > ${gwas_headers}.txt
ncol_id=`cat ${gwas_headers}.txt | awk '{if ($1 == "SNP" || $1 == "rs") print NR}'`
ncol_A1=`cat ${gwas_headers}.txt | awk '{if ($1 == "A1" || $1 == "allele1") print NR}'`
ncol_A2=`cat ${gwas_headers}.txt | awk '{if ($1 == "A2" || $1 == "allele0") print NR}'`
ncol_freq=`cat ${gwas_headers}.txt | awk '{if ($1 == "Freq" || $1 ~ "AF1") print NR}'`
ncol_b=`cat ${gwas_headers}.txt | awk '{if ($1 == "b" || $1 == "BETA") print NR}'`
ncol_se=`cat ${gwas_headers}.txt | awk '{if ($1 == "se" || $1 == "SE") print NR}'`
ncol_p=`cat ${gwas_headers}.txt | awk '{if ($1 == "p" || $1 == "P" || $1 == "p_lhr") print NR}'`
cat ${gwas} | awk -v OFS="\t" -v ncol_id=$ncol_id -v ncol_A1=$ncol_A1 -v ncol_A2=$ncol_A2 -v ncol_freq=$ncol_freq -v ncol_b=$ncol_b -v ncol_se=$ncol_se -v ncol_p=$ncol_p 'BEGIN{print "SNP\tA1\tA2\tfreq\tb\tse\tp\tn"}{if (NR >= 2) print $ncol_id, $ncol_A1, $ncol_A2, $ncol_freq, $ncol_b, $ncol_se, $ncol_p, "NA"}' > ${gwas_prefix}
rm -f ${gwas_headers}


#### Prepare eqtl file
echo -e "\n******** Prepare eqtl results ********\n"
eqtl_file_prefix="$(mktemp Temp/eqtl_file.XXXXXXXXXX)"
if [[ ${eqtl_format} == "MatrixEQTL" ]]
then
    # Make a BESD file from Matrix eQTL output
    ${smr} --eqtl-summary ${eqtl_file} --matrix-eqtl-format --thread-num ${threads} --make-besd --out ${eqtl_file_prefix}
elif [[ ${eqtl_format} == "tensorQTL" ]]
then
    #Make a BESD file from tensorQTL output
    if [[ ${eqtl_file} == *gz ]]
    then
        catcmd="zcat ${eqtl_file}"
    else
        catcmd="cat ${eqtl_file}"
    fi
    ${catcmd} | awk -v OFS="\t" '{if (NR >= 2) print $1, $2, $3, $7, $8}' | pigz -c -p ${threads} > ${eqtl_file_prefix}.gz
    ${smr} --eqtl-summary ${eqtl_file_prefix}.gz --fastqtl-nominal-format --thread-num ${threads} --make-besd --out ${eqtl_file_prefix}
    rm -f ${eqtl_file_prefix}.gz
elif [[ ${eqtl_format} == "fastqtl" ]]
then
    #Make a BESD file from fastqtl output
    ${smr} --eqtl-summary ${eqtl_file} --fastqtl-nominal-format --thread-num ${threads} --make-besd --out ${eqtl_file_prefix}
else
    echo "Error: We currently not support eQTL results from ${eqtl_format} !!!"
    exit_abnormal
    exit 1
fi

# Update coordinates of SNPs and genes, frequency of effect allele
echo -e "\n******** Update eqtl .esi .epi ********\n"
if [[ -f ${eqtl_map}.bim ]]
then
    myfreq="$(mktemp Temp/eqtl_myfreq.XXXXXXXXXX)"
    myesi="$(mktemp Temp/myesi.XXXXXXXXXX)"
    join -1 1 -2 2 <(zcat ${eqtl_file} | awk -v OFS="\t" '{if (NR >= 2) print $2, $1, $4, $8}' | sort -T Temp/ -k1,1) <(cat ${eqtl_map}.bim | sort -k2,2) | sed 's/ /\t/g' | awk -v OFS="\t" '{print $1,$8,$9,$3}' | sort -T Temp/ | uniq > ${myfreq}
    join -1 1 -2 2 <(zcat ${eqtl_file} | awk -v OFS="\t" '{if (NR >= 2) print $2, $1, $4, $8}' | sort -T Temp/ -k1,1) <(cat ${eqtl_map}.bim | sort -k2,2) | sed 's/ /\t/g' |  awk -v OFS="\t" '{print $5, $1,$6,$7,$8,$9,$4}' | sort -T Temp/ | uniq > ${myesi}
else
    echo "Error: Plink bim for eqtl mapping not exist!!!"
    exit_abnormal
    exit 1
fi
if [[ ${eqtl_bed} == *.gz ]]
then
    opencmd="zcat ${eqtl_bed}"
else
    opencmd="cat ${eqtl_bed}"
fi
myepi="$(mktemp Temp/myepi.XXXXXXXXXX)"
join -1 1 -2 4 <(zcat ${eqtl_file} | awk -v OFS="\t" '{if (NR >= 2) print $1, $2, $4, $8}' | sort -T Temp/ -k1,1) <(${opencmd} | awk -v OFS="\t" '{if (NR >= 2) print $1,$2,$3,$4}' | sort -T Temp/ -k4,4) | sed 's/ /\t/g' | awk '{if ($6 < $7) print $0"\t+"; else print $0"\t-"}' | awk -v OFS="\t" '{print $5, $1, "0", $6,$1, $8}' | sort -T Temp/ | uniq  >  ${myepi}
${smr} --beqtl-summary ${eqtl_file_prefix} --thread-num ${threads} --update-esi ${myesi}
${smr} --beqtl-summary ${eqtl_file_prefix} --thread-num ${threads} --update-epi ${myepi}
${smr} --beqtl-summary ${eqtl_file_prefix} --thread-num ${threads} --update-freq ${myfreq}
rm -f ${myesi}
rm -f ${myepi}
rm -f ${myfreq}

#### run SMR and HEIDI tests
echo -e "\n******** Run SMR test ********\n"
if [[ ${out} =~ "/" ]]
then
    outpath=${out%/*}
    mkdir -p ${outpath}
fi
${smr} --bfile ${bfile} --gwas-summary ${gwas_prefix} --beqtl-summary ${eqtl_file_prefix} --out ${out}.SMR --thread-num ${threads} --diff-freq 0.8 --diff-freq-prop 0.2 --heidi-mtd 1 --peqtl-smr 5e-6

sleep 30
if [[ -f ${out}.SMR.smr ]] 
then
    rm -rf ${myfreq}
    rm -rf ${myepi}
    rm -rf ${myesi}
    rm -rf ${bfile}
    rm -rf ${gwas_headers}
    rm -rf ${gwas_prefix}
    rm -rf ${eqtl_file_prefix}
fi
echo -e "Done!!\n"
