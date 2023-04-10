#!/usr/bin/env python3

'''
Take an annotation, breakpoint bedpe, and aligned bam
to produce a browser shot with annotations on top and
reads on the next panel that show the fusion

usage:
    python3 visualize_fusion.py \
        -o fusion_browser_shot.png \
        -a gencode.v38.annotation.gtf \
        -f prefix.breakpoints.groups.bed \
        -b prefix.mapped.bam
'''


import argparse
import os
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
import pickle


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--output', '-o', type=str,
        help='Output png filename'
    )
    parser.add_argument(
        '--annotation', '-a', type=str,
        help='Annotation gtf file'
    )
    parser.add_argument(
        '--fusion', '-f', type=str,
        help='A single line from the breakpoints groups file with the query fusion'
    )
    parser.add_argument(
        '--bam', '-b', type=str,
        help='Mapped bam file'
    )
    parser.add_argument(
        '--pickle', '-p', type=str, default='pbfusion_vis.pickle',
        help='Pickle filename [pbfusion_vis.pickle]'
    )
    return parser.parse_args()


def read_bedpe(bedpe):
    '''
    Takes the bedpe file and gives back:
        reads
        breakpoints
        gene ids
    '''
    with open(bedpe) as f:
        for line in f:
            if line.startswith('#'):
                continue
            sl = line.rstrip().split('\t')
            breakpoints = ((sl[0], sl[1]), (sl[3], sl[4]))
            reads = set(sl[11][6:].split(','))
            aux = sl[10].split(';')
            for field in aux:
                if field.startswith('GENE_IDS='):
                    gene_ids = field[9:].split(',')
    return reads, gene_ids, breakpoints


def read_gtf(gtf, gene_ids):
    gtf_dict = {}
    transcript_list = []
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            a = line.strip().split('\t')
            chromosome = a[0]
            type1 = a[2]
            if type1 in ['exon', 'CDS']:
                start = int(a[3])
                end = int(a[4])
                transcript = a[8].split(' transcript_id "')[1].split('"')[0]
                gene = a[8].split('gene_id "')[1].split('"')[0]
                if gene not in gene_ids:
                    continue
                if transcript not in gtf_dict:
                    gtf_dict[transcript] = []
                gtf_dict[transcript].append([chromosome, start, end, type1])

    for transcript, parts in gtf_dict.items():
        starts, ends, blockstarts, blockwidths, types = [], [], [], [], []
        for part in parts:
            starts.append(part[1])
            ends.append(part[2])
            blockstarts.append(part[1])
            blockwidths.append(part[2] - part[1])
            chromosome = part[0]
            types.append(part[3])
        transcript_list.append([
            chromosome, min(starts), max(ends),
            blockstarts, blockwidths, False, types])
    return transcript_list


def read_bam(bam, reads):
    pysam.set_verbosity(0)

    read_dict = {}
    with pysam.AlignmentFile(bam, 'rb', check_sq=False, require_index=False) as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if read.query_name in reads:
                if read.query_name not in read_dict:
                    read_dict[read.query_name] = []

                # make blockstarts and blockwidths
                blocks = []
                for b in read.get_blocks():
                    blocks.append((b[0], b[1]-b[0]))
                read_dict[read.query_name].append([read.reference_name, blocks])
    return read_dict


def plot_fusion(ref_genes, alignments, breakpoints, output):
    alignments = {k: v for (k, v) in alignments.items() if len(v) == 2}

    gene1_start = min([x[0][1][0][0] for x in list(alignments.values()) if x[0][0] == breakpoints[0][0]])
    gene2_end = max([x[2] for x in ref_genes if x[0] == breakpoints[1][0]])
    gene1_span = abs(int(breakpoints[0][1]) - gene1_start)
    gene2_span = abs(int(breakpoints[1][1]) - gene2_end)
    total_span = gene1_span + gene2_span
    gene1_prop = gene1_span / total_span
    gene2_prop = 1 - gene1_prop

    plt.figure(figsize=(8, 3))
    plt.style.use('clean')
    gene1 = plt.axes([0.1, 0.8, 0.8, 0.1])
    gene2 = plt.axes([0.1, 0.7, 0.8, 0.1])
    reads1 = plt.axes([0.1, 0.1, gene1_prop*0.8, 0.6])
    reads2 = plt.axes([0.1+(gene1_prop*0.8), 0.1, (gene2_prop*0.8), 0.6])

    orange, blue = (249/255, 157/255, 65/255), (19/255, 131/255, 198/255)

    # plot the gene annotations
    g1_y, g2_y = 0, 0
    for transcript in ref_genes:
        if transcript[0] == breakpoints[0][0]:
            panel = gene1
            g1_y += 1
            ypos = g1_y
            color = orange
        elif transcript[0] == breakpoints[1][0]:
            panel = gene2
            g2_y += 1
            ypos = g2_y
            color = blue
        start, stop = transcript[1], transcript[2]
        panel.plot([start, stop], [ypos]*2, lw=0.2, c=color, zorder=10)
        for i in range(len(transcript[3])):
            left, width, ftype = transcript[3][i], transcript[4][i], transcript[6][i]
            if ftype == 'exon':
                offset = 0.125
            else:
                offset = 0.25
            feature = Rect((left, ypos-offset), width, offset*2, lw=0, fc=color, zorder=12)
            panel.add_patch(feature)

    y = 0
    for read in list(alignments.values()):
        y += 1
        for locus in read:
            if locus[0] == breakpoints[0][0]:
                panel = reads1
                color = orange
                start, stop = locus[1][0][0], int(breakpoints[0][1])
            elif locus[0] == breakpoints[1][0]:
                panel = reads2
                color = blue
                start, stop = locus[1][0][0], locus[1][-1][0] + locus[1][-1][1]
            for block in locus[1]:
                left, width = block
                exon = Rect((left, y-0.25), width, 0.5, lw=0, fc=color, zorder=12)
                panel.add_patch(exon)
            panel.plot([start, stop], [y]*2, lw=0.2, c=color, zorder=10)

    gene1.set_ylabel('ASPSCR1', fontsize=5)
    gene2.set_ylabel('TFE3', fontsize=5)
    reads1.set_ylabel('Reads')

    gene1.set_xlim(gene1_start, int(breakpoints[0][1])+gene2_span)
    gene2.set_xlim(int(breakpoints[1][1])-gene1_span, gene2_end)
    reads1.set_ylim(0, len(alignments)+1)
    reads2.set_ylim(0, len(alignments)+1)
    reads1.set_xlim(gene1_start, int(breakpoints[0][1]))
    reads2.set_xlim(int(breakpoints[1][1]), gene2_end)

    reads1.spines['right'].set_visible(False)
    reads2.spines['left'].set_visible(False)
    for panel in [gene1, gene2, reads1, reads2]:
        panel.tick_params(
            bottom=False, labelbottom=False,
            left=False, labelleft=False,
            right=False, labelright=False,
            top=False, labeltop=False
        )

    if not output.endswith('.png'):
        output += '.png'
    plt.savefig(output, dpi=1000)


def main(args):
    reads, gene_ids, breakpoints = read_bedpe(args.fusion)
    if os.path.exists(args.pickle):
        with open(args.pickle, 'rb') as f:
            ref_genes = pickle.load(f)
            alignments = pickle.load(f)
    else:
        reads, gene_ids, breakpoints = read_bedpe(args.fusion)
        ref_genes = read_gtf(args.annotation, gene_ids)
        alignments = read_bam(args.bam, reads)
        with open(args.pickle, 'wb') as f:
            pickle.dump(ref_genes, f)
            pickle.dump(alignments, f)

    plot_fusion(ref_genes, alignments, breakpoints, args.output)

if __name__ == '__main__':
    args = parse_args()
    main(args)
