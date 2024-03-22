<h1 align="center"><img width="300px" src="pbfusion_logo.svg"/></h1>

<h1 align="center">pbfusion</h1>

Fusion gene caller for Iso-Seq sequencing data.

Authors: [Roger Volden](https://github.com/velociroger-pb), [Zev Kronenberg](https://github.com/zeeev), [Daniel Baker](https://github.com/dnbaker), [Khi Pin Chua](https://github.com/proteinosome)

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda) for information on Installation, Support, License, Copyright, and Disclaimer.


## Table of Contents ##
1. [Install](#install)
2. [Usage](#usage)
3. [Output](#output)
4. [Examples](#examples)
5. [Accessory scripts](#accessory)
6. [Help](#help)


## Install <a name="install"></a>
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pbfusion/README.html)

`pbfusion` can be installed from bioconda
```
conda install -c bioconda pbfusion
```
Binaries are also availible in the github releases.

## Usage <a name="usage"></a>

`pbfusion` has two subcommands: `pbfusion discover` and `pbfusion gff-cache`.

`pbfusion gff-cache` is not required, but recommended when running `pbfusion` multiple times. `pbfusion gff-cache` will serialize the input gtf/gff file and preprocess into exonic intervals ahead of time, which is fairly slow to do on the fly.
You can find GENCODE annotation files [here](https://www.gencodegenes.org/human/release_38.html).

```
Usage: pbfusion gff-cache [OPTIONS] --gtf <ReferenceAnnotation>

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
Identify fusion genes in aligned PacBio Iso-Seq data

Usage: pbfusion discover [OPTIONS] --gtf <REF> --output-prefix <OUTPUT> [FILE]...

Arguments:
  [FILE]...

Options:
  -b <ADDITIONAL_BAMS>
          Aligned Iso-Seq data in BAM format. Accepts a path to a bam, a url (if compiled with curl support), or a fofn (file-of-filenames) file with one filename or url per line
  -g, --gtf <REF>
          Reference gene annotations in GTF format. We also accept `gtf.bin` files as built by `pbfusion gff-cache`. This file must have `bin` as its suffix to be recognized. We also support gtf.bin.xz and gtf.bin.gz, compressed by xz and gzip, respectively. Recognition is based entirely on filename. Warning: the binary cached format has been altered since 0.3.3. You may need to re-generate your binary annotations.
  -o, --output-prefix <OUTPUT>
          Output prefix [default: none]
  -Q, --min-fusion-quality <MIN_FUSION_QUALITY>
          Determine the minimum fusion quality to emit. Choices: must be LOW or MEDIUM [default: MEDIUM]
  -t, --threads <THREADS>
          Number of threads. Defaults to available parallelism [default: 0]
  -c, --min-coverage <MIN_COVERAGE>
          Real-cell filtering for single-cell data. Iso-Seq reads annotated with zero "rc" tag value will be filtered. Assigns "low confidence" to fusion calls with read coverage below the minimum coverage threshold [default: 2]
  -i, --min-mean-identity <MIN_MEAN_IDENTITY>
          Assigns "low confidence" to fusion calls where the mean alignment identity is below the threshold [default: 0.85]
  -p, --min-mean-mapq <MIN_MEAN_MAPQ>
          Assigns "low confidence" to fusion calls where the mean mapq is below the threshold [default: 0.]
  -M, --min-fusion-read-fraction <MIN_FUSION_READ_FRACTION>
          Remove breakpoint pairs from groups if they have gene alignments which fewer than \[arg\] reads in group have [default: 0.25]
  -s, --max-variability <MAX_VARIABILITY>
          Assigns "low confidence" to fusion calls with the mean breakpoint distance is above the threshold [default: 1000]
  -a, --max-readthrough <MAX_READTHROUGH>
          Assigns "low confidence" to fusion calls spanning two genes below the readthrough threshold. [default: 100000]
  -m, --max-genes-in-event <MAX_GENES_IN_EVENT>
          Mark fusion groups involving > \[arg\] genes as low quality. This is a common source of false positives [default: 3]
  -r, --real-cell-filtering
      --allow-immune
          Permit fusion events identified involving primarily immunological genes and their pseudogenes. These are a common source of false positives and we mark them low-quality by default.
      --allow-mito
          Permit fusion events identified involving mitochondrial genes. These are a common source of false positives and we mark them low-quality by default.
      --prom-filter <PROM_FILTER>
          Filter rarer events involving genes with high numbers of fusion partners. These are a common source of false positives. Disable by setting `--prom-filter 0`. [default: 8]
  -v, --verbose...
          Enable verbose output
      --log-level <LOG_LEVEL>
          Alternative to repeated -v/--verbose: set log level via key.
          Values: "error", "warn" (default), "info", "debug", "trace".
          Enabling any level higher than "warn" also emits verbose output, including extra output files.
          If -v/--verbose is set, this option is ignored.
          Equivalence to -v/--verbose:
                => "WARN"
             -v => "INFO"
            -vv => "DEBUG"
           -vvv => "TRACE" [default: error]
  -h, --help
          Print help information (use `-h` for a summary)
  -V, --version
          Print version information

Copyright (C) 2004-2023     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.
```

## Output <a name="output"></a>

`pbfusion discover` produces **one output file designed for end users**: `{prefix}.breakpoints.groups.bed`. All other files are auxiliary and usually used for diagnostic purposes.

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
  - The cluster score is set to "MEDIUM" unless the group is considered low-quality, in which case it is "LOW"
    - Scores may be set to "LOW" depending on:
      - Low coverage
      - Low alignment identity
      - High breakpoint variability
      - Majority immunoglobulin genes
      - Mitochondrial genes
      - Having many gene partners across the dataset
  - You can `grep` the output for specific gene names or gene IDs if you're fishing for a specific fusion gene pair


### Filtering options

Filtering is now primarily done through adjusting the allowed score [default "MEDIUM"].
The primary filtering options to reduce false positives are `--min-coverage` (`-c`), `--max-readthrough` (`-a`), `--min-mean-identity` (`-i`), `--min-mean-mapq` (`-p`), `--min-fusion-read-fraction` (`-M`), and `--max-variability` (`-s`).
`--min-coverage` just filters out breakpoints based on their read support, where the default value is 2 to filter out singletons.
`--max-readthrough` is used to discard reads that align to two genes next to each other in the genome [default 100kb].
`--min-mean-identity` will assign low confidence to fusions with mean mapping identity lower than the threshold [0.85]. Sometimes in hard to map regions of the genome, the aligner will force an alignment through a region and incur a high edit distance.
`--min-mean-mapq` default is set to `0.` because even with high alignment identity, some genes have close relatives where the aligner will assign a MAPQ of 0.
`--min-fusion-read-fraction` is used to filter long chains of genes.
An example of this would be an IG alignment, where maybe 100 reads align to various annotated regions (eg. IGHA1, IGHV3-23, IGHV3-7, IGHG3, IGHJ5, IGHJ4, IGHJ3, IGHGP, IGHJ2, IGHJ6).
Given a read coverage of these genes like [100, 100, 90, 90, 20, 10, 10, 5, 5, 5], we would by default filter genes with coverage lower than 25% of the total read count for this fusion.
With that filter, you're left with this coverage: [100, 100, 90, 90], which with the read filtering brings it down to [80, 80, 70, 70].
Because this fusion still has >3 genes in it, it would get filtered out.
`--max-variability` allows you to filter based on breakpoint variability [default 1000].

As of v0.4.0, the default behavior is to mark entries with simple majority of immune genes as `LOW`.
We use the GENCODE `gene_type` field to classify annotations as immune.
This filter can be disabled by setting the `--allow-immune` option.
Additionally, we mark entries with mitochondrial genes as `LOW`, which can be disabled with the `--allow-mito` option.
Lastly, we have implemented a promiscuity filter to help decrease false positives.
This filter works by taking all of the `MEDIUM` fusion entries and tracking how many different gene partners each gene has.
If `Gene_A` has entries with genes `B`, `C`, ..., `K`, then `Gene_A` will be subject to the promiscuity filter [default is 8 gene partners].
For genes with many partners, we calculate the expected coverage for these entries as `sum(read_coverage) // n_partners`.
Entries that do not pass this expected coverage will get marked as `LOW`.


## Examples <a name="examples"></a>

Serializing the input gtf file:
```
pbfusion gff-cache \
    --gtf gencode.v38.annotation.gtf \
    --gtf-out gencode.v38.annotation.gtf.bin
```

Running `pbfusion` on aligned reads:
```
pbfusion discover \
    --reference-gtf gencode.v38.annotation.gtf.bin \
    --output-prefix isoseq \
    --threads 8 \
    isoseq.aligned.bam
```

You may find substantial space savings by compressing your annotation bin file.
pbfusion will accept the smaller file. Being substantially smaller (5M vs 140M), this makes artefact management much easier.

```
> xz -c -9 --extreme gencode.v38.annotation.gtf.bin > gencode.v38.annotation.gtf.bin.xz
> ls -oh gencode.v38.annotation.gtf.bin* | awk '{print $4, $NF}'
140M gencode.v38.annotation.gtf.bin
5.6M gencode.v38.annotation.gtf.bin.xz
```


## Accessory scripts <a name="accessory"></a>
There are two scripts included in this repo for convenience: `visualize_fusion.py` and `extract_tag.py`.  


### Visualize fusion
`visualize_fusion.py` will produce a genome browser plot when given an annotation, a single-line BEDPE file, and a mapped BAM.
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


### Extract tag
**This script is deprecated now that cell barcodes are emitted automatically**
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


### Changelog
Changelog - PacBio Fusion Detection - pbfusion

## v0.4.1 3/22/24
### Changes
- Bug fix: consistent ordering between gene ids and gene names when dealing with ambiguity.
- Bug fix: rely on breakpoint ordering in visualization script.
- Improvement: gene names are assigned automatically in visualization script.
- Improvement: visualization script can now take SAM files as input.
- Improvement: visualization script can now take a breakpoint number (these aren't always unique, so it's best to still give a single entry).

## v0.4.0 11/28/23
### Changes
- Improvement: Improved ambiguous gene resolution.
- Improvement: read lzma-compressed annotations in text and binary format.
- Improvement: Support multiple bam/fofn (File-Of-FileNames) inputs.
- Improvement: reduced false-positive fusions.
- Improvement: Additional metadata for putative fusion candidates.
- Improvement: Updated cached binary format. WARNING: this is a breaking change. Old binary annotation files will need to be re-generated.
- Improvement: Maintain read ordering in reporting fusions.
- Bug fix: Update reported breakpoints so that the interval described is the last-in-read and the second is first-soft-clipped.
- Improvement: Add gene names to annotation in order they occur in reads.
- Improvement: Report all breakpoint coordinates for reads in a breakpoint cluster, add as a new tag.
- Improvement: Update read-through annotation such that we define it based on the breakpoint distance instead of the distance between genes.
- Alteration: reduced chaining distance for clustering breakpoints.
- Alteration: emit sense-antisense events by default.

## v0.3.1 8/17/23
### Changes
- Bug fix: resolving query position sorting.

## v0.3.0 7/25/23
### Changes
- Bug fix: resolving strand-related breakpoint errors.
- Bug fix: GTF parsing of overlapping annotations.
- Bug fix: annoting gtf records without parent gene entries.

## v0.2.3 7/11/23
### Changes
- Add sample name, sequencing platform, and read group information to output BED file header.

## v0.2.2 6/29/23
### Changes
- Command-line interface, filtering, and formatting changes.
- Option to filter by average mapping quality on either side of a breakpoint.
- Option to filter reads from groups based on gene coverage (as a %-age of total reads for that group).
- Filters out entries where the number of genes is 1 or >3.
- We now mark events by quality, and then filter by quality levels instead of directly filtering events.
  These events can be emitted with `--min-fusion-quality LOW` instead of the default `MEDIUM`.
- Simplified output. Fewer categories, simpler identifiers, and a new header format that is easier to parse.
- Number of reads filter now checks the `im` bam tag for clustered input. Now, the number of original reads is used for filtering instead of the number of clusters.

## v0.2.1: 6/8/23
### Changes
- Add options to filter segments by length and error rate.

## v0.2.0: 5/22/23
### Changes
- Join pbfusion and pbgffcache into one executable `pbfusion` with subcommands `pbfusion discover` and `pbfusion gff-cache` corresponding to the previous.

## v0.1.2: 5/19/23
### Changes
- Add edit distance per segment as a feature of Interval for filtering.

## v0.1.1: 4/26/23
### Changes
- Add cell barcodes if found, and original read names if reads are clustered to output file. The first tag set in `CB`, `XC`, and `CR` is used.
- Add `--max-mean-breakpoint-distance` flag, which marks breakpoint groups as low quality if the mean distance of breakpoints in a cluster exceeds this threshold.

## v0.1.0
### Changes
- Initial release

## Help <a name="help"></a>

Please direct support/help/bug questions to [Github issue](https://github.com/PacificBiosciences/pbbioconda/issues)

## Disclaimer

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.


