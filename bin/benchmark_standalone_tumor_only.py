#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from pathlib import Path
import csv
import os

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


def filter_germline(all_called_variants):
    """
    Filter germline variants based on VAF (tumor_only mode only)
    """

    dict = {}

    for key, val in all_called_variants.items():
        try:
            df = val.loc[val['VAF'].astype(float) < 0.3]
            dict[key] = df
        except ValueError:
            continue
    return dict


def load_somatic(vcf_path):
    """
    Create somatic dataframes from BAMSurgeon VCF files (spike-in variants)
    """
    vcf_files = [path for path in Path(vcf_path).glob('*.vcf.gz')]
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_somatic = pd.DataFrame()

    for file in vcf_files:
        samplename = file.stem.split('.')[0]

        if 'bamsurgeon' in samplename:
            samplename = samplename.replace('bamsurgeon_', '')
        if 'spiked_snv' in samplename:
            samplename = samplename.replace('_spiked_snv', '')
        if 'spiked_indel' in samplename:
            samplename = samplename.replace('_spiked_indel', '')

        for variant in VCF(file):
            tmp = pd.DataFrame(data=
            [
                [samplename, variant.CHROM,
                 variant.POS, variant.REF,
                 variant.ALT[0], variant.INFO['VAF']
                 ]
            ], columns=df_cols
            )
            df_somatic = pd.concat([df_somatic, tmp])
    df_somatic.set_index('sample')
    return df_somatic


def load_callers(vcf_path):
    """
    Create dictionary with variants from variant callers VCF files
    """

    vcf_files = [path for path in Path(vcf_path).glob('**/*.vcf.gz')]

    variants = {}
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_vardict_snv = pd.DataFrame(columns=df_cols)
    df_vardict_indel = pd.DataFrame(columns=df_cols)
    df_vardict = pd.DataFrame(columns=df_cols)

    df_mutect_snv = pd.DataFrame(columns=df_cols)
    df_mutect_indel = pd.DataFrame(columns=df_cols)
    df_mutect = pd.DataFrame(columns=df_cols)

    df_varscan_snv = pd.DataFrame(columns=df_cols)
    df_varscan_indel = pd.DataFrame(columns=df_cols)
    df_varscan = pd.DataFrame(columns=df_cols)

    df_freebayes_snv = pd.DataFrame(columns=df_cols)
    df_freebayes_indel = pd.DataFrame(columns=df_cols)
    df_freebayes = pd.DataFrame(columns=df_cols)

    df_lofreq_snv = pd.DataFrame(columns=df_cols)
    df_lofreq_indel = pd.DataFrame(columns=df_cols)
    df_lofreq = pd.DataFrame(columns=df_cols)

    for file in vcf_files:

        name = file.name

        if 'vardictjava' in name:
            samplename = file.stem.replace('vardictjava_', '').split('.')[0]
            for variant in VCF(file):
                if variant.is_snp:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_vardict_snv = pd.concat([df_vardict_snv, tmp])
                elif variant.is_indel:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_vardict_indel = pd.concat([df_vardict_indel, tmp])


        elif 'gatk4_mutect2' in name:
            samplename = file.stem.replace('gatk4_mutect2_', '').split('.')[0]
            for variant in VCF(file):
                if variant.is_snp:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('AF')[0][0]
                         ]
                    ], columns=df_cols
                    )
                    df_mutect_snv = pd.concat([df_mutect_snv, tmp])
                elif variant.is_indel:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.format('AF')[0][0]
                         ]
                    ], columns=df_cols
                    )
                    df_mutect_indel = pd.concat([df_mutect_indel, tmp])


        elif 'varscan2' in name:
            samplename = file.stem.replace('varscan2_', '').split('.')[0]
            for variant in VCF(file):
                tmp = pd.DataFrame(data=
                [
                    [samplename, variant.CHROM,
                     variant.POS, variant.REF,
                     variant.ALT[0], float(variant.format('FREQ')[0].strip('%')) / 100
                     ]
                ], columns=df_cols
                )
                df_varscan_snv = pd.concat([df_varscan_snv, tmp])


        elif 'freebayes' in name:
            samplename = file.stem.replace('freebayes_', '').split('.')[0]
            for variant in VCF(file):
                if variant.is_snp:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_freebayes_snv = pd.concat([df_freebayes_snv, tmp])
                elif variant.is_indel:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_freebayes_indel = pd.concat([df_freebayes_indel, tmp])


        elif 'lofreq' in name:
            samplename = file.stem.replace('lofreq_', '').split('.')[0]
            for variant in VCF(file):
                if variant.is_snp:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_lofreq_snv = pd.concat([df_lofreq_snv, tmp])
                if variant.is_indel:
                    tmp = pd.DataFrame(data=
                    [
                        [samplename, variant.CHROM,
                         variant.POS, variant.REF,
                         variant.ALT[0], variant.INFO['AF']
                         ]
                    ], columns=df_cols
                    )
                    df_lofreq_indel = pd.concat([df_lofreq_indel, tmp])


    df_vardict_snv.set_index('sample')
    df_vardict_indel.set_index('sample')
    df_vardict = pd.concat([df_vardict_snv, df_vardict_indel])

    df_mutect_snv.set_index('sample')
    df_mutect_indel.set_index('sample')
    df_mutect = pd.concat([df_mutect_snv, df_mutect_indel])

    df_varscan_snv.set_index('sample')
    df_varscan_indel.set_index('sample')
    df_varscan = pd.concat([df_varscan_snv, df_varscan_indel])

    df_freebayes_snv.set_index('sample')
    df_freebayes_indel.set_index('sample')
    df_freebayes = pd.concat([df_freebayes_snv, df_freebayes_indel])

    df_lofreq_snv.set_index('sample')
    df_lofreq_indel.set_index('sample')
    df_lofreq = pd.concat([df_lofreq_snv, df_lofreq_indel])

    variants = {
        'vardict_snvs': df_vardict_snv,
        'vardict_indels': df_vardict_indel,
        'vardict_all': df_vardict,

        'mutect_snvs': df_mutect_snv,
        'mutect_indels': df_mutect_indel,
        'mutect_all': df_mutect,

        'varscan_snvs': df_varscan_snv,
        'varscan_indels': df_varscan_indel,
        'varscan_all': df_varscan,

        'freebayes_snvs': df_freebayes_snv,
        'freebayes_indels': df_freebayes_indel,
        'freebayes_all': df_freebayes,

        'lofreq_snvs': df_lofreq_snv,
        'lofreq_indels': df_lofreq_indel,
        'lofreq_all': df_lofreq
    }

    return variants


def calculate_spikein(called, spiked_snv,
                      spiked_indel=pd.DataFrame(
                          columns=['sample', 'chrom', 'pos', 'REF', 'ALT']
                      )
                      ):
    """
    Calculate the number of spiked-in variants called in the input variant caller
    """
    dict = {}

    for key, val in called.items():
        if 'snv' in key:
            df = val.merge(spiked_snv, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key] = df
        if 'indel' in key:
            df = val.merge(spiked_indel, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key] = df
        if 'all' in key:
            df = val.merge(pd.concat([spiked_snv, spiked_indel]),
                           on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key] = df
    # if germinal is not None:
    # df_spiked_germinal = df.merge(df_neat, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")

    return dict


def calculate_performance(true_called, all_called, spiked_snv, spiked_indel=pd.DataFrame()):
    """
    Calculate performance metrics for input the variant caller
    """
    cols = ['TP', 'FP', 'FN', 'TPR', 'PPV', 'FDR']
    dict = {}
    for key, val in true_called.items():
        tp = len(val)
        fp = len(all_called[key]) - tp
        if 'snv' in key:
            fn = len(spiked_snv) - tp
        if 'indel' in key:
            fn = len(spiked_indel) - tp
        if 'all' in key:
            fn = len(spiked_indel) + len(spiked_snv) - tp
        try:
            tmp = pd.DataFrame(data=
            [
                [tp, fp, fn,
                 tp / (tp + fn), tp / (tp + fp), fp / (fp + tp)]
            ], columns=cols
            )
            dict[key] = tmp
        except ZeroDivisionError:
            continue
    return dict


def plot_performance(performance, output, type):
    """
    Plot performance for the variant callers
    """

    # Set Black and white cmap
    from matplotlib import lines, markers

    # line cyclers adapted to colourblind people
    from cycler import cycler
    line_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                   cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
    marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                     cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                     cycler(marker=["^", "P", ".", "1", "+", "x", "."]))

    if type == 'snv' or type == 'indel':
        # Create figure object and store it in a variable called 'fig'
        fig = plt.figure(figsize=(5, 5))

        # Add axes object to our figure that takes up entire figure
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_prop_cycle(marker_cycler)

        if type == 'snv':
            for key, df in performance.items():
                if 'snv' in key:
                    ax.plot(df['PPV'], df['TPR'], markersize=15, label=key.split('_')[0])
        elif type == 'indel':
            for key, df in performance.items():
                if 'indel' in key:
                    ax.plot(df['PPV'], df['TPR'], markersize=15, label=key.split('_')[0])

        # Add the x and y-axis labels
        ax.set_xlabel(r'PPV', labelpad=10, fontsize=20)
        ax.set_ylabel(r'Sensitivity', labelpad=10, fontsize=20)

        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)

        # Add legend to plot
        ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False, fontsize=13)

        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.grid(axis='both', linestyle='--')

        plt.savefig(output + '/plots/benchmark_' + type + '.png', dpi=350, transparent=False,
                    bbox_inches='tight')

    elif type == 'both':
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
        for key, df in performance.items():
            if 'snv' in key:
                ax1.plot(df['PPV'], df['TPR'], markersize=15, label=key.split('_')[0])

        for key, df in performance.items():
            if 'indel' in key:
                ax2.plot(df['PPV'], df['TPR'], markersize=15, label=key.split('_')[0])

        for key, df in performance.items():
            if 'all' in key:
                ax3.plot(df['PPV'], df['TPR'], markersize=15, label=key.split('_')[0])

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
        plt.savefig(output + '/plots/benchmark_' + type + '.png', dpi=350, transparent=False,
                    bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser(description='Generate plots for artificial mutation benchmark')

    parser.add_argument('-n', '--normal', help='Normal samples')

    parser.add_argument('-t', '--mut_type', default='snv', help='Type of mutation spiked-in the samples.'
                                                                'Can be: snv, indel, both')

    parser.add_argument('-s', '--snv',
                        help='Folder with VCF files with spiked in somatic snv')

    parser.add_argument('-v', '--variant_calling',
                        help='Folder with VCF files with variant callers output')
    parser.add_argument('-i', '--indel',
                        help='Folder with VCF files with spiked in somatic indel')
    parser.add_argument('-o', '--output',
                        help='Output folder')

    args = parser.parse_args()

    if not os.path.exists(args.output + "/spiked_variants"):
        os.makedirs(args.output + "/spiked_variants")
        os.makedirs(args.output + "/all_called_variants")
        os.makedirs(args.output + "/true_called_variants")
        os.makedirs(args.output + "/plots")

    if args.normal:
        # Load pseudo-germinal variants (generated from NEAT)
        df_normal = load_germinal(args.normal)
        df_normal.to_csv("germinal_variants.txt", sep="\t")

    # Load ground-truth somatic variants (spiked-in from BAMSurgeon)
    if args.mut_type == 'snv':
        spiked_snv = load_somatic(args.snv)
        spiked_snv.to_csv(args.output + "/spiked_variants/somatic_spiked_snvs.txt", sep="\t")
        spiked_snv.to_excel(args.output + "/spiked_variants/somatic_spiked_snvs.xlsx")
    elif args.mut_type == 'indel':
        spiked_indel = load_somatic(args.indel)
        spiked_indel.to_csv(args.output + "/spiked_variants/somatic_spiked_indels.txt", sep="\t")
        spiked_indel.to_excel(args.output + "/spiked_variants/somatic_spiked_indels.xlsx")
    else:
        spiked_snv = load_somatic(args.snv)
        spiked_snv.to_csv(args.output + "/spiked_variants/somatic_spiked_snvs.txt", sep="\t")
        spiked_snv.to_excel(args.output + "/spiked_variants/somatic_spiked_snvs.xlsx")

        spiked_indel = load_somatic(args.indel)
        spiked_indel.to_csv(args.output + "/spiked_variants/somatic_spiked_indels.txt", sep="\t")
        spiked_indel.to_excel(args.output + "/spiked_variants/somatic_spiked_indels.xlsx")

    # Load variants from variant callers vcfs
    all_called_variants = load_callers(args.variant_calling)

    # Filter germline variants in tumor_only mode
    all_called_variants_f = filter_germline(all_called_variants)

    # Save variants for each caller in txt files
    for key, val in all_called_variants_f.items():
        val.to_csv(args.output + "/all_called_variants/" + key + ".txt", sep="\t")
        val.to_excel(args.output + "/all_called_variants/" + key + ".xlsx")

    # Calculate variants spiked for each caller
    if args.mut_type == 'snv':
        true_called_variants = calculate_spikein(all_called_variants, spiked_snv)
    else:
        true_called_variants = calculate_spikein(all_called_variants, spiked_snv, spiked_indel)

    for key, val in true_called_variants.items():
        val.to_csv(args.output + "/true_called_variants/" + key + ".txt", sep="\t")
        val.to_excel(args.output + "/true_called_variants/" + key + ".xlsx")

    # Calculate performance for each caller
    if args.mut_type == 'snv':
        performance = calculate_performance(true_called_variants, all_called_variants, spiked_snv)
    else:
        performance = calculate_performance(true_called_variants, all_called_variants, spiked_snv, spiked_indel)

    p_ = pd.DataFrame()
    for key, val in performance.items():
        val['Caller'] = key
        val = val.set_index('Caller')
        p_ = pd.concat([p_, val])
    p_.to_csv(args.output + "/performance" + ".txt", sep="\t")
    p_.to_excel(args.output + "/performance" + ".xlsx")

    # Plot performance
    plot_performance(performance, args.output, args.mut_type)


if __name__ == '__main__':
    main()
