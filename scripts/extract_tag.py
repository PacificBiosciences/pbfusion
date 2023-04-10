#!/usr/bin/env python3

'''
Usage:
    python3 extract_tag.py \
        -b reads.bam \
        -t CB

Optional:
    -i include.txt (list of read names)

Emits a tab delimited file with the read names and the
specified tag
'''

import argparse
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bam', '-b', required=True, type=str,
        help='Input bam file'
    )
    parser.add_argument(
        '--tag', '-t', required=True, type=str, default='CB',
        help='Tag that you want extracted [CB]'
    )
    parser.add_argument(
        '--include', '-i', type=str,
        help='List of read names to include (one per line)'
    )
    return parser.parse_args()

def read_include(include):
    to_include = set()
    with open(include) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            to_include.add(line)
    return to_include

def extract_tags(bam, tag, include=set()):
    pysam.set_verbosity(0)
    filter_reads = bool(include)
    with pysam.AlignmentFile(bam, 'rb', check_sq=False, require_index=False) as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if filter_reads and read.query_name not in include:
                continue
            try:
                extracted_tag = read.get_tag(tag)
            except KeyError:
                print(f"Couldn't find tag: {tag}", file=sys.stderr)
                exit(1)
            print(f'{read.query_name}\t{extracted_tag}')

def main(args):
    if args.include:
        include = read_include(args.include)
        extract_tags(args.bam, args.tag, include)
    else:
        extract_tags(args.bam, args.tag)

if __name__ == '__main__':
    args = parse_args()
    main(args)
