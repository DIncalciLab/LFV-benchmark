#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import matplotlib.font_manager as fm
from cyvcf2 import VCF
from sklearn import metrics


def load_germinal(vcf):
    """
    Create dataframes from NEAT VCF files (pseudo-germinal variants)
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT"]
    df_neat = pd.DataFrame()

    for file in vcf:
        print(file)
        samplename = file.replace('neat_', '').replace('.vcf.gz', '')
        for variant in VCF(file):
            df_neat = df_neat.append(
                [
                    [samplename, variant.CHROM,
                     variant.POS, variant.REF,
                     variant.ALT[0]
                     ]
                ]
            )

    df_neat.columns = df_cols
    return df_neat


def load_ground_truth(vcf):
    """
    Create grount truth dataframes from BAMSurgeon VCF files (spike-in variants)
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_bamsurgeon = pd.DataFrame()

    for file in vcf:
        samplename = file.replace('bamsurgeon_', '').replace('.vcf.gz', '')
        for variant in VCF(file):
            df_bamsurgeon = df_bamsurgeon.append(
                [
                    [samplename, variant.CHROM,
                     variant.POS, variant.REF,
                     variant.ALT[0], variant.INFO['VAF']
                     ]
                ]
            )
    df_bamsurgeon.columns = df_cols
    return df_bamsurgeon


def load_vardict(vcf):
    """
    Create dataframes from VarDict VCF files
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_vardict_snv = pd.DataFrame()
    df_vardict_indel = pd.DataFrame()

    for file in vcf:
        samplename = file.replace('vardictjava_', '').replace('.vcf.gz', '')
        for variant in VCF(file):
            if variant.is_snp:
                df_vardict_snv = df_vardict_snv.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ]
                )
            elif variant.is_indel:
                df_vardict_indel = df_vardict_indel.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ]
                )

    df_vardict_snv.columns = df_cols
    df_vardict_indel.columns = df_cols

    return df_vardict_snv, df_vardict_indel


def load_mutect(vcf):
    """
    Create dataframes from Mutect2 VCF files
    """

    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_mutect_snv = pd.DataFrame()
    df_mutect_indel = pd.DataFrame()

    for file in vcf:
        samplename = file.replace('mutect_', '').replace('.vcf.gz', '')
        for variant in VCF(file):
            if variant.is_snp:
                df_mutect_snv = df_mutect_snv.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('AF')[0][0]
                         ]
                    ]
                )
            elif variant.is_indel:
                df_mutect_indel = df_mutect_indel.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('AF')[0][0]
                         ]
                    ]
                )

    df_mutect_snv.columns = df_cols
    df_mutect_indel.columns = df_cols

    return df_mutect_snv, df_mutect_indel


def load_varscan(vcf):
    """
    Create dataframes from VarScan2 VCF files
    """

    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_varscan_snv = pd.DataFrame()
    df_varscan_indel = pd.DataFrame()

    for file in vcf:
        samplename = file.replace('varscan_', '').replace('.vcf.gz', '')
        for variant in VCF(file):
            if variant.is_snp:
                df_varscan_snv = df_varscan_snv.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('FREQ')[0]
                         ]
                    ]
                )
            elif variant.is_indel:
                df_varscan_indel = df_varscan_indel.append(
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('FREQ')[0]
                         ]
                    ]
                )

    df_varscan_snv.columns = df_cols
    df_varscan_indel.columns = df_cols

    return df_varscan_snv, df_varscan_indel


def calculate_spikein(df, df_bamsurgeon, df_neat):
    """
    Calculate the number of spiked-in variants called in the input variant caller
    """

    df_spiked = df.merge(df_bamsurgeon, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
    df_spiked_germinal = df.merge(df_neat, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")

    return df_spiked, df_spiked_germinal


def calculate_performance(df, df_spiked, df_spiked_germinal, df_truth):
    """
    Calculate performance metrics for input the variant caller
    """

    tp = len(df_spiked)
    fp = len(df) - len(df_spiked_germinal) - tp
    fn = len(df_truth) - tp

    performance = pd.DataFrame()

    performance = performance.append([
        [tp, fp, fn, tp / (tp + fn), tp / (tp + fp), fp / (fp + tp)]
    ])

    performance.columns = ['TP', 'FP', 'FN', 'TPR', 'PPV', 'FDR']

    return performance


def plot_performance(vardict, mutect, varscan):
    """
    Plot performance for the input variant caller
    """

    # Set Black and white cmap
    from matplotlib import lines, markers

    d = {'VarDict': vardict,
         'Mutect2': mutect,
         'VarScan2': varscan}

    # line cyclers adapted to colourblind people
    from cycler import cycler
    line_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                   cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
    marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                     cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                     cycler(marker=["^", "P", ".", "1", "+", "x", "."]))

    # Create figure object and store it in a variable called 'fig'
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 5))

    # Edit the major and minor ticks of the x and y axes
    # ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    # ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    # ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')

    ax1.set_prop_cycle(marker_cycler)
    ax2.set_prop_cycle(marker_cycler)
    ax3.set_prop_cycle(marker_cycler)

    # Plot the sample
    for key, df in d.items():
        ax1.plot(df['PPV'], df['TPR'], markersize=15, label=key)

    # for key, df in d_ind.items():
    #    ax2.plot(df['PPV'], df['TPR'], markersize=15, label=key)

    # for key, df in d_tot.items():
    #    ax3.plot(df['PPV'], df['TPR'], markersize=15, label=key)

    # Edit the major and minor tick locations of x and y axes
    # ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    # ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

    # ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

    # Add the x and y-axis labels
    ax1.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    ax1.set_ylabel(r'Sensitivity', labelpad=10, fontsize=20)

    ax2.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    # ax2.set_ylabel(r'Sensitivity', labelpad=10, fontsize=15)

    ax3.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    # ax3.set_ylabel(r'Sensitivity', labelpad=10, fontsize=15)

    # Set the axis <limits
    ax1.set_xlim(-0.1, 1.1)
    ax1.set_ylim(-0.1, 1.1)

    ax2.set_xlim(-0.1, 1.1)
    ax2.set_ylim(-0.1, 1.1)

    ax3.set_xlim(-0.1, 1.1)
    ax3.set_ylim(-0.1, 1.1)

    ax1.text(-0.2, 1.2, 'a)', fontsize=20)
    ax2.text(-0.2, 1.2, 'b)', fontsize=20)
    ax3.text(-0.2, 1.2, 'c)', fontsize=20)

    # Hide the top and right spines of the axis
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)

    ax1.grid(axis='both', linestyle='--')
    ax2.grid(axis='both', linestyle='--')
    ax3.grid(axis='both', linestyle='--')

    # Add legend to plot
    # ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=10)
    # ax2.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=10)
    ax3.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=15)

    # Add titles to subplot
    ax1.set_title('Performance for SNVs', fontsize=20)
    ax2.set_title('Performance for INDELs', fontsize=20)
    ax3.set_title('Performance for SNVs and INDELs', fontsize=20)

    # Save figure
    plt.savefig('benchmark.png', dpi=350, transparent=False, bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser(description='Generate plots for artificial mutation benchmark')

    parser.add_argument('-n', '--normal', nargs='+', help='Normal samples')

    parser.add_argument('-b', '--somatic', nargs='+', help='Somatic samples')

    parser.add_argument('-v', '--vardict', nargs='+',
                        help='VarDict calls')
    parser.add_argument('-m', '--mutect', nargs='+',
                        help='Mutect2 calls')
    parser.add_argument('-s', '--varscan', nargs='+',
                        help='VarScan2 calls')
    parser.add_argument('-f', '--freebayes', nargs='+',
                        help='FreeBayes calls')
    parser.add_argument('-k', '--strelka', required=True, nargs='+',
                        help='Strelka2 calls')
    parser.add_argument('-lofreq', '--lofreq', required=True, nargs='+',
                        help='LoFreq calls')

    args = parser.parse_args()

    if (args.normal):
        # Load pseudo-germinal variants (generated from NEAT)
        df_normal = load_germinal(args.normal)
        df_normal.to_csv("germinal_variants.txt", sep="\t")

    # Load ground-truth somatic variants (spiked-in from BAMSurgeon)
    df_somatic = load_ground_truth(args.somatic)
    df_somatic.to_csv("somatic_spikein_variants.txt", sep="\t")

    # Load VarDict variants
    df_vardict_snv, df_vardict_indel = load_vardict(args.vardict)
    df_vardict_tot = df_vardict_snv.append(df_vardict_indel)
    df_vardict_tot.to_csv("vardict_variants.txt", sep="\t")

    # Load Mutect2 variants
    df_mutect_snv, df_mutect_indel = load_mutect(args.mutect)
    df_mutect_tot = df_mutect_snv.append(df_mutect_indel)
    df_mutect_tot.to_csv("mutect_variants.txt", sep="\t")

    # Load VarScan2 variants
    df_varscan_snv, df_varscan_indel = load_varscan(args.varscan)
    df_varscan_tot = df_varscan_snv.append(df_varscan_indel)
    df_varscan_tot.to_csv("varscan_variants.txt", sep="\t")

    # Calculate spiked-in variants for each caller
    df_spikein_vardict, df_spikein_vardict_germinal = calculate_spikein(df_vardict_tot, df_somatic, df_normal)
    df_spikein_mutect, df_spikein_mutect_germinal = calculate_spikein(df_mutect_tot, df_somatic, df_normal)
    df_spikein_varscan, df_spikein_varscan_germinal = calculate_spikein(df_varscan_tot, df_somatic, df_normal)

    # Calculate performance for each caller
    df_performance_vardict = calculate_performance(df_vardict_tot, df_spikein_vardict, df_spikein_vardict_germinal,
                                                   df_normal)
    df_performance_vardict.to_csv("vardict_performance.txt", sep="\t")

    df_performance_mutect = calculate_performance(df_mutect_tot, df_spikein_mutect, df_spikein_mutect_germinal, df_normal)
    df_performance_mutect.to_csv("mutect_performance.txt", sep="\t")

    df_performance_varscan = calculate_performance(df_varscan_tot, df_spikein_varscan, df_spikein_varscan_germinal,
                                                   df_normal)
    df_performance_varscan.to_csv("varscan_performance.txt", sep="\t")

    # Plot performance
    plot_performance(df_performance_vardict, df_performance_mutect, df_performance_varscan)


if __name__ == '__main__':
    main()
