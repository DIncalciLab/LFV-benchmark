#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from cyvcf2 import VCF


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


def load_somatic(vcf):
    """
    Create somatic dataframes from BAMSurgeon VCF files (spike-in variants)
    """
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_somatic = pd.DataFrame()

    for file in vcf:
        if 'bamsurgeon' in file:
            samplename = file.replace('bamsurgeon_', '').replace('.vcf.gz', '')
        else:
            samplename = os.path.basename(vcf).split('.')[0]
        for variant in VCF(file):
            df_somatic = df_somatic.append(
                [
                    [samplename, variant.CHROM,
                     variant.POS, variant.REF,
                     variant.ALT[0], variant.INFO['VAF']
                     ]
                ]
            )
    df_somatic.columns = df_cols
    return df_somatic


def load_callers(vcf):
    """
    Create dictionary with variants from variant callers VCF files
    """

    variants = {}
    df_cols = ["sample", "chrom", "pos", "REF", "ALT", "VAF"]

    df_vardict_snv = pd.DataFrame()
    df_vardict_indel = pd.DataFrame()
    df_vardict = pd.DataFrame()

    df_mutect_snv = pd.DataFrame()
    df_mutect_indel = pd.DataFrame()
    df_mutect = pd.DataFrame()

    df_varscan_snv = pd.DataFrame()
    df_varscan_indel = pd.DataFrame()
    df_varscan = pd.DataFrame()

    df_freebayes_snv = pd.DataFrame()
    df_freebayes_indel = pd.DataFrame()
    df_freebayes = pd.DataFrame()

    df_strelka_snv = pd.DataFrame()
    df_strelka_indel = pd.DataFrame()
    df_strelka = pd.DataFrame()

    df_lofreq_snv = pd.DataFrame()
    df_lofreq_indel = pd.DataFrame()
    df_lofreq = pd.DataFrame()

    for file in vcf:

        if 'vardictjava' in file:
            samplename = file.replace('vardictjava_', '').replace('.vcf', '')
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
            df_vardict = pd.concat([df_vardict_snv, df_vardict_indel])

        elif 'gatk_mutect2' in file:
            samplename = file.replace('gatk_mutect2_', '').replace('.vcf', '')
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
            df_mutect = pd.concat([df_mutect_snv, df_mutect_indel])

        elif 'varscan2' in file:
            samplename = file.replace('varscan2_', '').replace('.vcf', '')
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
            df_varscan = pd.concat([df_varscan_snv, df_varscan_indel])

        elif 'freebayes' in file:
            samplename = file.replace('freebayes_', '').replace('.vcf.gz', '')
            for variant in VCF(file):
                if variant.is_snp:
                    df_freebayes_snv = df_freebayes_snv.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], variant.INFO['AF']
                             ]
                        ]
                    )
                elif variant.is_indel:
                    df_freebayes_indel = df_freebayes_indel.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], variant.INFO['AF']
                             ]
                        ]
                    )

            df_freebayes_snv.columns = df_cols
            df_freebayes_indel.columns = df_cols
            df_freebayes = pd.concat([df_freebayes_snv, df_freebayes_indel])


        elif 'strelka' in file:
            if 'snvs' in file:
                samplename = file.replace('strelka_somatic_', '').replace('.somatic_snvs.vcf.gz', '')
                for variant in VCF(file):
                    refCounts = variant.REF + 'U'
                    altCounts = variant.ALT[0] + 'U'

                    tier1RefCounts = variant.format(refCounts)[0][0]
                    tier1AltCounts = variant.format(altCounts)[1][0]

                    vaf = (tier1AltCounts) / (tier1AltCounts + tier1RefCounts)

                    df_strelka_snv = df_strelka_snv.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], vaf
                             ]
                        ]
                    )
            elif 'indels' in file:
                samplename = file.replace('strelka_somatic_', '').replace('.somatic_indels.vcf.gz', '')
                for variant in VCF(file):
                    tier1RefCounts = variant.format('TAR')[0]
                    tier1AltCounts = variant.format('TIR')[0]

                    vaf = (tier1AltCounts) / (tier1AltCounts + tier1RefCounts)
                    df_strelka_indel = df_strelka_indel.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], vaf
                             ]
                        ]
                    )

            df_strelka_snv.columns = df_cols
            df_strelka_indel.columns = df_cols
            df_strelka = pd.concat([df_strelka_snv, df_strelka_indel])

        elif 'lofreq' in file:
            if 'snvs' in file:
                samplename = file.replace('strelka_somatic_', '').replace('.somatic_snvs.vcf.gz', '')
                for variant in VCF(file):
                    refCounts = variant.REF + 'U'
                    altCounts = variant.ALT[0] + 'U'

                    tier1RefCounts = variant.format(refCounts)[0][0]
                    tier1AltCounts = variant.format(altCounts)[1][0]

                    vaf = (tier1AltCounts) / (tier1AltCounts + tier1RefCounts)

                    df_lofreq_snv = df_lofreq_snv.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], vaf
                             ]
                        ]
                    )
            elif 'indels' in file:
                samplename = file.replace('strelka_somatic_', '').replace('.somatic_indels.vcf.gz', '')
                for variant in VCF(file):
                    tier1RefCounts = variant.format('TAR')[0]
                    tier1AltCounts = variant.format('TIR')[0]

                    vaf = (tier1AltCounts) / (tier1AltCounts + tier1RefCounts)
                    df_lofreq_indel = df_lofreq_indel.append(
                        [
                            [samplename, variant.CHROM,
                             variant.POS, variant.REF,
                             variant.ALT[0], vaf
                             ]
                        ]
                    )

            df_lofreq_snv.columns = df_cols
            df_lofreq_indel.columns = df_cols
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

        'strelka_snvs': df_strelka_snv,
        'strelka_indels': df_strelka_indel,
        'strelka_all': df_strelka,

        'lofreq_snvs': df_lofreq_snv,
        'lofreq_indels': df_lofreq_indel,
        'lofreq_all': df_lofreq
    }

    return variants


def calculate_spikein(called, somatic, germinal=None):
    """
    Calculate the number of spiked-in variants called in the input variant caller
    """
    dict = {}

    for key, val in called.items():
        if 'all' in key:
            all = called.merge(somatic, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key].append(all)
        elif 'snvs' in key:
            snvs = called.merge(somatic, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key].append(snvs)
        elif 'indels' in key:
            indels = called.merge(somatic, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")
            dict[key].append(indels)
    # if germinal is not None:
    # df_spiked_germinal = df.merge(df_neat, on=['sample', 'chrom', 'pos', 'REF', 'ALT'], how="inner")

    return dict


def calculate_performance(spiked, somatic, germinal = None):
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


def plot_performance(spike_variants, somatic, germinal = None):
    """
    Plot performance for the variant callers
    """

    # Set Black and white cmap
    from matplotlib import lines, markers

    #d = {'VarDict': vardict,
         #'Mutect2': mutect,
         #'VarScan2': varscan}

    # line cyclers adapted to colourblind people
    from cycler import cycler
    line_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                   cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
    marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                     cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                     cycler(marker=["^", "P", ".", "1", "+", "x", "."]))

    # Create figure object and store it in a variable called 'fig'
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, figsize=(20, 5))

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

    parser.add_argument('-s', '--somatic', nargs='+', help='Somatic samples')

    parser.add_argument('-v', '--variant_calling', nargs='+',
                        help='VCF files with variant callers output')

    args = parser.parse_args()

    if args.normal:
        # Load pseudo-germinal variants (generated from NEAT)
        df_normal = load_germinal(args.normal)
        df_normal.to_csv("germinal_variants.txt", sep="\t")

    # Load ground-truth somatic variants (spiked-in from BAMSurgeon)
    df_somatic = load_somatic(args.somatic)
    df_somatic.to_csv("somatic_spike_variants.txt", sep="\t")

    # Load variants from variant callers vcfs
    dict_variants = load_callers(args.variant_calling)

    # Save variants for each caller in txt files
    for key, val in dict_variants.items():
        val.to_csv(key + ".txt", sep="\t")

    # Calculate variants spiked for each caller
    spike_variants = calculate_spikein(dict_variants)

    # Calculate spiked-in variants for each caller
    # df_spikein_vardict, df_spikein_vardict_germinal = calculate_spikein(df_vardict_tot, df_somatic, df_normal)
    # df_spikein_mutect, df_spikein_mutect_germinal = calculate_spikein(df_mutect_tot, df_somatic, df_normal)
    # df_spikein_varscan, df_spikein_varscan_germinal = calculate_spikein(df_varscan_tot, df_somatic, df_normal)

    # Calculate performance for each caller
    performance = calculate_performance(spike_variants, df_somatic)
    """
    df_performance_vardict = calculate_performance(df_vardict_tot, df_spikein_vardict, df_spikein_vardict_germinal,
                                                   df_normal)
    df_performance_vardict.to_csv("vardict_performance.txt", sep="\t")

    df_performance_mutect = calculate_performance(df_mutect_tot, df_spikein_mutect, df_spikein_mutect_germinal,
                                                  df_normal)
    df_performance_mutect.to_csv("mutect_performance.txt", sep="\t")

    df_performance_varscan = calculate_performance(df_varscan_tot, df_spikein_varscan, df_spikein_varscan_germinal,
                                                   df_normal)
    df_performance_varscan.to_csv("varscan_performance.txt", sep="\t")
    """
    # Plot performance
    plot_performance(df_performance_vardict, df_performance_mutect, df_performance_varscan)


if __name__ == '__main__':
    main()
