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
        samplename = file.name.split('.')[0]
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



def load_ground_truth(vcf_snv, vcf_indel):
    """
    Create grount truth dataframes from BAMSurgeon VCF files (spike-in variants)
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT"]

    df_groundtruth_snv = pd.DataFrame()
    df_groundtruth_indel = pd.DataFrame()


    for file in vcf_snv:
        samplename = file.name.split('.')[0]
        for variant in VCF(file):
            df_groundtruth_snv = df_groundtruth_snv.append(
                [
                    [samplename, variant.CHROM,
                    variant.POS, variant.REF,
                    variant.ALT[0], variant.INFO['VAF']
                    ]
                ]
            )
            
    df_groundtruth_snv.columns = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    for file in vcf_indel:
        samplename = file.name.split('.')[0]
        for variant in VCF(file):
            df_groundtruth_indel = df_groundtruth_indel.append(
                [
                    [samplename, variant.CHROM,
                    variant.POS, variant.REF,
                    variant.ALT[0], variant.INFO['VAF']
                    ]
                ]
            )
            
    df_groundtruth_indel.columns = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    return df_groundtruth_snv, df_groundtruth_indel

def load_vardict(vcf):
    """
    Create dataframes from VarDict VCF files
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_vardict_snv = pd.DataFrame()
    df_vardict_indel = pd.DataFrame()

    for file in vardict_files:
        samplename = file.name.split('.')[0]
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
        samplename = file.name.split('.')[0]
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

    for file in varscan_files:
        samplename = file.name.split('.')[0]
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


def calculate_spikein(df, df_truth, df_germinal):
    """
    Calculate the number of spiked-in variants called in the input variant caller
    """

    df_spiked = df.merge(df_truth, on=['sample','chrom', 'pos', 'REF', 'ALT'], how="inner")
    df_spiked_germinal = df.merge(df_germinal, on=['sample','chrom', 'pos', 'REF', 'ALT'], how="inner")

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
        [tp, fp, fn, tp/(tp+fn), tp/(tp+fp), fp/(fp+tp)]
    ])

    performance.columns = ['TP', 'FP', 'FN', 'TPR', 'PPV', 'FDR']

    performance.to_csv(output)
    return performance

def plot_performance(vardict, mutect, varscan):
    """
    Plot performance for the input variant caller
    """

    #Set Black and white cmap
    from matplotlib import lines, markers

    d = {'VarDict': vardict, 
         'Mutect2': mutect, 
         'VarScan2': varscan}

    # line cyclers adapted to colourblind people
    from cycler import cycler
    line_cycler   = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                    cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
    marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                    cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                    cycler(marker=["^", "P", ".", "1", "+", "x", "."]))

    
     # Create figure object and store it in a variable called 'fig'
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,5))


    # Edit the major and minor ticks of the x and y axes
    #ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    #ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='out')
    #ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')

    ax1.set_prop_cycle(marker_cycler)
    ax2.set_prop_cycle(marker_cycler)
    ax3.set_prop_cycle(marker_cycler)

    # Plot the sample
    for key, df in d_snv.items():
        ax1.plot(df['PPV'], df['TPR'], markersize=15, label=key)
        
    for key, df in d_ind.items():
        ax2.plot(df['PPV'], df['TPR'], markersize=15, label=key)
        
    for key, df in d_tot.items():
        ax3.plot(df['PPV'], df['TPR'], markersize=15, label=key)

    # Edit the major and minor tick locations of x and y axes
    #ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    #ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

    #ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))

    # Add the x and y-axis labels
    ax1.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    ax1.set_ylabel(r'Sensitivity', labelpad=10, fontsize=20)

    ax2.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    #ax2.set_ylabel(r'Sensitivity', labelpad=10, fontsize=15)

    ax3.set_xlabel(r'PPV', labelpad=10, fontsize=20)
    #ax3.set_ylabel(r'Sensitivity', labelpad=10, fontsize=15)

    # Set the axis <limits
    ax1.set_xlim(-0.1, 1.1)
    ax1.set_ylim(-0.1, 1.1)

    ax2.set_xlim(-0.1, 1.1)
    ax2.set_ylim(-0.1, 1.1)

    ax3.set_xlim(-0.1, 1.1)
    ax3.set_ylim(-0.1, 1.1)

    ax1.text(-0.2,1.2,'a)', fontsize=20)
    ax2.text(-0.2,1.2,'b)', fontsize=20)
    ax3.text(-0.2,1.2,'c)', fontsize=20)

    # Hide the top and right spines of the axis
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)

    ax1.grid(axis = 'both', linestyle='--')
    ax2.grid(axis = 'both', linestyle='--')
    ax3.grid(axis = 'both', linestyle='--')


    # Add legend to plot
    #ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=10)
    #ax2.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=10)
    ax3.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False, fontsize=15)

    #Add titles to subplot
    ax1.set_title('Performance for SNVs', fontsize=20)
    ax2.set_title('Performance for INDELs', fontsize=20)
    ax3.set_title('Performance for SNVs and INDELs', fontsize=20)


    # Save figure
    plt.savefig('benchmark.png', dpi=350, transparent=False, bbox_inches='tight')


    def main():

        parser = argparse.ArgumentParser(description='Generate plots for artificial mutation benchmark')

        #parser.add_argument('-n', '--neat',       required=True, help='Pseudo-germinal variants generated from NEAT')

        parser.add_argument('-g', '--bamsurgeon', required=True, help='Mutations spiked-in from BamSurgeon')

        #parser.add_argument('-v', '--vardict',    required=True, help='VarDict VCF with spiked-in artificial mutations')
        #parser.add_argument('-s', '--varscan',    required=True, help='VarScan2 VCF with spiked-in artificial mutations')
        #parser.add_argument('-m', '--mutect',     required=True, help='Mutect2 VCF with spiked-in artificial mutations')

        args = parser.parse_args()

        #Load pseudo-germinal variants (generated from NEAT)
        #df_germinal = load_germinal(args.neat)

        df_germinal.read_excel("test.xlsx")
        
        #Load ground-truth variants (spiked-in from BAMSurgeon)
        df_truth = load_ground_truth(args.bamsurgeon)

        #Load VarDict variants

        #df_vardict = load_vardict(vcf_vardict)
        #Load Mutect2 variants
        #df_mutect = load_mutect(vcf_mutect)

        #Load VarScan2 variants
        #df_varscan = load_varscan(vcf_varscan)


        #Calculate spiked-in variants for each caller
        #df_spikein_vardict, df_spikein_vardict_germinal = calculate_spikein(df_vardict, df_truth, df_germinal)
        #df_spikein_mutect, df_spikein_mutect_germinal = calculate_spikein(df_mutect, df_truth, df_germinal)
        #df_spikein_varscan, df_spikein_varscan_germinal = calculate_spikein(df_varscan, df_truth, df_germinal)

        #Calculate performance for each caller
        #df_performance_vardict = calculate_performance(df_vardict, df_spikein_vardict, df_spikein_vardict_germinal, df_truth)
        #df_performance_mutect = calculate_performance(df_mutect, df_spikein_mutect, df_spikein_mutect_germinal, df_truth)
        #df_performance_varscan = calculate_performance(df_varscan, df_spikein_varscan, df_spikein_varscan_germinal, df_truth)

        #Plot performance
        #plot_performance(df_performance_vardict, df_performance_mutect, df_performance_varscan)

    if __name__ == '__main__':
        main()