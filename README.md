<h1 align="left"><img width="300px" src="logo.png"/></h1>

# pbfusion

Fusion gene detection tool for PacBio Iso-Seq data generated using HiFi sequencing. Can be used on bulk and single-cell Iso-Seq data.

authors: Roger Volden, Zev Kronenberg, Daniel Baker, Khi Pin Chua


## Table of Contents ##
1. [Install](#install)
2. [Usage](#usage)
3. [Output](#output)
4. [Examples](#examples)
5. [Accessory scripts](#accessory)


## Install <a name="install"></a>

`pbfusion` can be installed from bioconda
```
conda install -c bioconda pbfusion
```


## Usage <a name="usage"></a>

`pbfusion` has two primary executables: `pbfusion` and `gffcache`.

`gffcache` is not required, but recommended when running `pbfusion` multiple times. `gffcache` will serialize the input GTF/GFF file and preprocess into exonic intervals ahead of time, which is fairly slow to do on the fly.

```
Usage: gffcache [OPTIONS] --gtf <ReferenceAnnotation>

Options:
  -g, --gtf <ReferenceAnnotation>            Input GTF file
  -b, --gtf-out <BinaryReferenceAnnotation>  Output binary GTF file [default: *]
  -v, --verbose...                           Enable verbose output
  -h, --help                                 Print help information
  -V, --version                              Print version information
```

`pbfusion` requires two input files and one option:

1. Aligned Iso-Seq HiFi data in a BAM file format. The data should be aligned with `pbmm2` using the `ISOSEQ` preset and the `--sort` flag enabled. `pbfusion` accepts Iso-Seq reads or polished transcripts (the output of `isoseq3 cluster`).
2. Reference gene annotations in GTF format. This GTF file *must* match the reference genome used for alignments.
3. The output prefix. `pbfusion` writes multiple files, prefixed with the user specified string.

```
A tool for detecting fusion genes in aligned PacBio RNA data

Usage: pbfusion [OPTIONS] --bam <FILE> --reference-gtf <REF> --output-prefix <OUTPUT>

Options:
  -b, --bam <FILE>
          Aligned FLNC or transcript BAM file
  -g, --reference-gtf <REF>
          Input annotations in gtf format
  -o, --output-prefix <OUTPUT>
          Output bed file prefix [default: none]
  -d, --breakpoint-group-distance <MAX_CLUST_DIST>
          Maximum allowed distance (in bp) for clustering breakpoints [default: 500]
  -t, --threads <THREADS>
          Number of threads [default: 1]
  -m, --min-coverage <MIN_COVERAGE>
          Filter out breakpoint groups with coverage lower than this number of reads (exclusive) [default: 2]
  -v, --verbose...
          Verbose output
          Apply once for info-level messages and additional output information
          Apply twice or more for debug logs
          If --log-level is set, this option takes precedence.
      --log-level <LOG_LEVEL>
          Alternative to repeated -v/--verbose: set log level via key.
          Values: "error", "warn" (default), "info", "debug", "trace".
          Enabling any level higher than "warn" also emits verbose output, including extra output files.
          If -v/--verbose is set, this option is ignored.
          Equivalence to -v/--verbose:
                => "warn"
             -v => "info"
            -vv => "debug"
           -vvv => "trace" [default: error]
  -T, --fusion-readthrough-threshold <fusion_readthrough_threshold>
          Distance to use to determine if a transcript aligned to two nearby genes on the same strand
          is a read-through event or a canonical fusion. [default: 10000]
  -F, --filter-groups-failing-any
          Filter out groups with breakpoint pairs failing any category.
          Standard behavior is to filter out groups with pairs failing all categories, which is more permissive.
  -R, --emit-readthrough
          To emit ReadThrough events. These are not emitted by default
  -S, --emit-sense-antisense
          To emit SenseAntisense events. These are not emitted by default
  -A, --emit-unannotated-exons
          To emit unannotated exon events. These are not emitted by default
  -O, --emit-overlapping
          To emit overlapping gene events. These are not emitted by default
  -a, --emit-all-breakpoints
          Emit all breakpoints. Enables --emit-unannotated-exons, --emit-sense-antisense, --emit-overlapping, and --emit-readthrough
  -n, --novel-exon-discovery
          Whether to discover "novel exons"
  -h, --help
          Print help information (use `--help` for more detail)
  -V, --version
          Print version information

Copyright (C) 2004-2023     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.
```

By default, readthrough and sense-antisense chimeras, and overlapping gene events, and events involving unannotated exons are filtered out.

They can be emitted by enabling relevant `--emit-{category}` flags, or all can be enabled with `-a/--emit-all-breakpoints`.


## Output <a name="output"></a>

`pbfusion` produces **one output file designed for end users**: `{prefix}.breakpoints.groups.bed`. All other files are auxiliary and usually used for diagnostic purposes.

If verbose output is enabled (`-v`), five output files sharing the same prefix are emitted.


| File                                 | Description |
| ------------------------------------ | ----------- |
| {prefix}.breakpoints.bed             | All detected breakpoints, BED format, aux output         |
| {prefix}.transcripts                 | All transcripts with a breakpoint, plain text, aux output          |
| {prefix}.breakpoints.groups.bed      | Clustered breakpoint calls, BEDPE format, main output            |
| {prefix}.unannotated.bed             | unannotated aligned segments, BED format, aux output         |
| {prefix}.unannotated.clusters.bed    | clustered unannotated aligned segments, BED format, aux output      |


### Clustered breakpoint calls file format

The clustered breakpoint call file is [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) file formatted including header lines.

Notes:
  - The cluster score is set to a constant "MEDIUM"
  - You can `grep` the output for specific gene names or gene IDs if you're fishing for a specific fusion gene pair


### Filtering options

The primary filtering options to reduce false positives are `--min-coverage` (`-m`) and `--fusion-readthrough-threshold` (`-T`).
`--min-coverage` just filters out breakpoints based on their read support, where the default value is 2 to filter out singletons.
`--fusion-readthrough-threshold` is used to discard reads that align to two genes next to each other in the genome.
By default, this is set to a modest 10kb, but for more stringent fusion gene calling, we recommend setting this option to 100kb.
Another option to use for more stringent fusion gene calling is `--filter-groups-failing-any`.
This means that breakpoint pairs that fail any of our internal checks (readthrough, strand switch, unannotated exons, or overlapping genes), that pair will not get emitted.


## Examples <a name="examples"></a>

Serializing the input gtf file:
```
gffcache \
    --gtf gencode.v38.annotation.gtf \
    --gtf-out gencode.v38.annotation.gtf.bin
```

Running `pbfusion` on aligned reads:
```
pbfusion \
    --bam isoseq.aligned.bam \
    --reference-gtf gencode.v38.annotation.gtf.bin \
    --output-prefix isoseq \
    --threads 8
```


## Accessory scripts <a name="accessory"></a>
There are two scripts included in this repo for convenience: `visualize_fusion.py` and `extract_tag.py`.  


### Visualize fusion
`visualize_fusion.py` will produce a genome browser shot when given an annotation, a single-line BEDPE file, and a mapped BAM.
Optionally, you can specify a path to a pickle file for easier rerunning when refining a figure.

Usage:
```
python3 visualize_fusion.py \
    -o fusion_browser_shot.png \
    -a gencode.v38.annotation.gtf \
    -f isoseq.breakpoints.groups.bed \
    -b isoseq.mapped.bam
```

Visualization dependencies:
 - matplotlib
 - pysam

The read names for the y-axis labels are done manually, so they will need to be changed in the source code when run.


### Extract tag
The main use case for this script is to output associated cell barcodes for fusion gene reads.
`extract_tag.py` will take a BAM file, BAM tag, and a list of read names (can be taken from the output BEDPE and edited for one readname per line).
The output is tab delimited with the read name and its associated cell barcode.
This can be used to extract any BAM tag, but it will look for the `CB` tag by default.

Usage:
```
python3 extract_tag.py \
    -b isoseq.aligned.bam \
    -t CB \
    -i readnames.txt \
    >cell_associations.tsv
```

This script will be deprecated once cell barcodes are emitted automatically.
