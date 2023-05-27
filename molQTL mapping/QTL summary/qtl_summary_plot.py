import numpy as np
import pandas as pd
import os
from pathlib import Path
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import hsv_to_rgb
import seaborn as sns
import scipy.stats as stats
import subprocess
import qtl.plot
import qtl.annotation
import qtl.norm
import qtl.genotype as gt
import argparse

parser = argparse.ArgumentParser(description='Plotting QTL summary information')
parser.add_argument('phenotype', help='')
parser.add_argument('working_path', help='')
parser.add_argument('color_settings', help='')
args = parser.parse_args()



### define functions
def plot_qtl_summary(qtl_summary_df, subtype='',prefix='e', normalized=False, ylim=None):
    colors = colors_df.loc[qtl_summary_df.index, 'color_hex']
    ax = qtl.plot.setup_figure(2,2,xspace=[1,0.25])
    if not normalized:
        ax.scatter(qtl_summary_df['samples'],
                   qtl_summary_df['egenes']/qtl_summary_df['detected_genes'],
                   c=colors, edgecolor='k', s=30, label='Expression', lw=0.66, clip_on=False)
    else:
        ax.scatter(qtl_summary_df['samples'],
                   qtl_summary_df['egenes'+subtype]/qtl_summary_df['detected_genes'+subtype],
                   c=colors, edgecolor='k', s=30, label='Expression', lw=0.66, clip_on=False)
    qtl.plot.format_plot(ax, fontsize=12)
    if ylim is None:
        ax.set_ylim([0, ax.get_ylim()[1]])
    else:
        ax.set_ylim(ylim)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(min_n_ticks=4, nbins=4))
    ax.set_xlim([0, ax.get_xlim()[1]])
    ax.spines['bottom'].set_position(('outward', 4))
    ax.spines['left'].set_position(('outward', 4))
    ax.set_xlabel('Samples', fontsize=14)
    ax.set_ylabel('{}Genes/tested_Genes'.format(prefix), fontsize=14)
    return ax

def plot_ah_summary(summary_df, prefix='e', width=0.9):
    ratio_df = (summary_df[summary_df.columns[:4]].T / summary_df['{}Genes'.format(prefix)]).T
    ax = qtl.plot.setup_figure(4,2,xspace=[0.75,0.75])
    c = summary_df.columns[:4][::-1]
    ratio_df[c].plot(kind='bar', stacked=True, ax=ax, width=width, color=sns.color_palette("Blues_r", 4))
    mean_color = hsv_to_rgb([0.02,0.8,0.8])
    ax2 = ax.twinx()
    x = np.arange(summary_df.shape[0])
    ax2.plot(x, summary_df['mean'], '.', c=mean_color, label='Mean')
    ax2.spines['right'].set_color(mean_color)
    ax2.tick_params(axis='y', colors=mean_color)
    qtl.plot.format_plot(ax, fontsize=12)
    qtl.plot.format_plot(ax2, fontsize=12, hide=['left', 'top', 'bottom'])
    ax2.set_ylim([0.5,2.5])
    ax.set_ylabel('Proportion of {}Genes'.format(prefix), fontsize=14)
    ax2.set_ylabel('Mean {}QTLs/gene'.format(prefix), fontsize=14, labelpad=8)#, c=mean_color)
    ylim = [0,14000]
    ylim = [0,1]
    d = ylim[1]-ylim[0]
    s = 1.06
    b = ylim[1] - s*d
    ax.set_ylim(ylim)
    ax.tick_params(axis='x', which='major', pad=10, labelsize=10)
    ax.scatter(x, b*np.ones(summary_df.shape[0]),s=20,c=colors_df.loc[summary_df.index, 'color_hex'], clip_on=False)
    ax.spines['left'].set_position(('outward', 4))
    ax2.spines['right'].set_position(('outward', 4))
    ax.spines['bottom'].set_bounds([-width/2, summary_df.shape[0]-1+width/2])
    ax.set_xlim([-0.66, summary_df.shape[0]-1+0.66])
    h1,_ = ax.get_legend_handles_labels()
    h2,_ = ax2.get_legend_handles_labels()
    h = h2 + h1
    labels = ['1', '2', '3', '≥ 4', 'Mean']
    h = h1
    labels = ['1', '2', '3', '≥ 4']
    ax.legend(h[::-1], labels, loc='upper left', title='{}QTLs/gene'.format(prefix),labelspacing=0.33, handletextpad=0.5)
    ax.set_xticklabels(colors_df.loc[indep_eqtl_summary_df.index].index,rotation=-90, ha='center', fontsize=6)


### Input color settings
colors_df = pd.read_csv(args.color_settings, sep="\t", index_col=0)
qtl_summary_file=args.working_path+"/"+args.phenotype+"_qtl_summary.txt"
eqtl_sum_pdf=args.working_path+"/"+args.phenotype+"_qtl_summary.pdf"
qtl_summary_df = pd.read_csv(qtl_summary_file, sep='\t', index_col=0, dtype=str)
qtl_summary_df=qtl_summary_df.apply(pd.to_numeric)
### Plot sample size vs eGenes
if args.phenotype=="expression":
    prefix="e"
elif args.phenotype=="splicing":
    prefix="s"
elif args.phenotype=="lncRNA":
    prefix="lnc"
elif args.phenotype=="exon":
    prefix="ex"
elif args.phenotype=="APA":
    prefix="3a"
ax = plot_qtl_summary(qtl_summary_df, prefix=prefix)
plt.savefig(eqtl_sum_pdf)

indep_summary_file=args.working_path+"/"+args.phenotype+"_indep_summary.txt"
indep_eqtl_summary_df = pd.read_csv(indep_summary_file, sep='\t', index_col=0, dtype=str)
indep_eqtl_summary_df=indep_eqtl_summary_df.apply(pd.to_numeric)
### Plot indep eQTLs
allele_het_pdf=args.working_path+"/"+args.phenotype+"_indep_summary.pdf"
plot_ah_summary(indep_eqtl_summary_df, prefix=prefix)
plt.savefig(allele_het_pdf,bbox_inches='tight')