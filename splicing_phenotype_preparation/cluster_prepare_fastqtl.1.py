from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import gzip
import contextlib
from datetime import datetime
import tempfile
import shutil
import glob
from sklearn.decomposition import PCA

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def get_ann_bed(annotation_bed, feature='gene'):
    """
    Parse genes from GTF, create placeholder DataFrame for BED output
    """
    chrom = []
    start = []
    end = []
    gene_id = []
    with open(annotation_bed, 'r') as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            chrom.append(row[0])

            # TSS: gene start (0-based coordinates for BED)
            if row[3]=='+':
                start.append(np.int64(row[1])-1)
                end.append(np.int64(row[1]))
            elif row[3]=='-':
                start.append(np.int64(row[2])-1)  # last base of gene
                end.append(np.int64(row[2]))
            else:
                raise ValueError('Strand not specified.')

            gene_id.append(row[4])

    bed_df = pd.DataFrame(
        data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id},
        columns=['chr', 'start', 'end', 'gene_id'],
        index=gene_id)
    return bed_df


def write_bed(bed_df, output_name):
    """Write DataFrame to BED"""
    bgzip = subprocess.Popen('bgzip -c > '+output_name,
        stdin=subprocess.PIPE, shell=True)
    bed_df.to_csv(bgzip.stdin, sep='\t', index=False)
    stdout, stderr = bgzip.communicate()
    subprocess.check_call('tabix -f '+output_name, shell=True)


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run leafcutter clustering, prepare for FastQTL')
    parser.add_argument('junc_files_list', help='File with paths to ${sample_id}.junc files')
    parser.add_argument('exons', help='Exon definitions file, with columns: chr, start, end, strand, gene_id, gene_name')
    parser.add_argument('genes_bed', help='Collapsed gene annotation in bed format')
    parser.add_argument('prefix', help='Prefix for output files (sample set ID)')
    parser.add_argument('--min_clu_reads', default='50', type=str, help='Minimum number of reads supporting each cluster')
    parser.add_argument('--min_clu_ratio', default='0.001', type=str, help='Minimum fraction of reads in a cluster that support a junction')
    parser.add_argument('--max_intron_len', default='500000', type=str, help='Maximum intron length')
    parser.add_argument('--leafcutter_dir', default='/opt/leafcutter',
                        help="leafcutter directory, containing 'clustering' directory")
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('-chr', '--chr_list', default='.', help='Chromosome list (only taking autosome into account)')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] leafcutter clustering')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    print('  * decompressing and renaming junc files')
    with open(args.junc_files_list) as f:
        junc_files = f.read().strip().split('\n')

    junc_dir = os.path.join(args.output_dir)
    if not os.path.exists(junc_dir):
        os.mkdir(junc_dir)

    # select autosomes
    chrs = pd.read_csv(args.chr_list, header=None)
    chrs = [str(x) for x in chrs[0].tolist()]

    sample_ids = []
    for f in junc_files:
        sample_id = os.path.split(f)[1].split('_rg')[0]
        sample_ids.append(sample_id)
        tmp_jun_file = pd.read_csv(f, sep = "\t", header=None)
        tmp_jun_file.rename(columns = {0:'chr'},inplace = True)
        tmp_jun_file_chr=tmp_jun_file[tmp_jun_file['chr'].isin(chrs)]
        #tmp_jun_file["chr"] = pd.to_numeric(tmp_jun_file["chr"])]
        #tmp_jun_file_chr=tmp_jun_file[tmp_jun_file['chr'].apply(lambda x: isinstance(x, (int, np.int64)))]
        out_jun_file=junc_dir+'/'+sample_id+'.junc'
        tmp_jun_file_chr.to_csv(out_jun_file, index=False, header = False, sep='\t')
        #shutil.copy2(f, os.path.join(junc_dir, sample_id)+'.junc')

    #subprocess.check_call('gunzip -f '+os.path.join(junc_dir, '*.junc.gz'), shell=True)
    junc_files = sorted([os.path.join(junc_dir, i+'.junc') for i in sample_ids])

    print('  * running leafcutter clustering')
    # generates ${prefix}_perind_numers.counts.gz and ${prefix}_perind.counts.gz
    with tempfile.NamedTemporaryFile(dir=args.output_dir) as tmp:
        with open(tmp.name, 'w') as f:
            f.write('\n'.join(junc_files)+'\n')
        subprocess.check_call(
            'python2 '+os.path.join(args.leafcutter_dir, 'clustering', 'leafcutter_cluster.py' \
                +' --juncfiles '+tmp.name \
                +' --outprefix '+args.output_dir + '/' +args.prefix \
                +' --minclureads '+args.min_clu_reads \
                +' --mincluratio '+args.min_clu_ratio \
                +' -s True --maxintronlen '+args.max_intron_len), shell=True)

    print('  * compressing outputs')
    subprocess.check_call('gzip -f '+os.path.join(args.output_dir,args.prefix+'_pooled'), shell=True)
    subprocess.check_call('gzip -f '+os.path.join(args.output_dir,args.prefix+'_refined'), shell=True)

    print('  * mapping clusters to genes')
    subprocess.check_call(
        'Rscript' \
            +' map_clusters_to_genes.R' \
            +' '+os.path.join(args.output_dir, args.prefix+'_perind.counts.gz') \
            +' '+args.exons \
            +' '+args.output_dir+'/'+args.prefix + '.leafcutter.clusters_to_genes.txt', shell=True)
print('done')
