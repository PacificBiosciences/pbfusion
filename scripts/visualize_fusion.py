#!/usr/bin/env python3

"""
Take an annotation, breakpoint bedpe, and aligned bam
to produce a browser shot with annotations on top and
reads on the next panel that show the fusion

usage:
    python3 visualize_fusion.py \
        -o fusion_browser_shot.png \
        -a gencode.v38.annotation.gtf \
        -f prefix.breakpoints.groups.bed \
        -b prefix.mapped.bam
"""


import argparse
from dataclasses import dataclass
from enum import Enum
import logging
from pathlib import Path
import pickle
import sys
from typing import Dict, List, NamedTuple, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
import pysam


logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="[%(asctime)s] %(levelname)s - %(message)s")
logger = logging.getLogger(Path(__file__).stem)

plt.style.use(Path(__file__).parent / "clean.mplstyle")

# The colors to use for each side of the fusion
GENE1_COLOR = "#DF1995"  # PacBio Bright Magenta
GENE2_COLOR = "#1383C6"  # PacBio Bright Blue


# margins to leave open around the panels
PLOT_MARGIN = 0.05
# minimum and maximum proportion of the y plot area to give to the reference transcripts
MIN_REF_TRANSCRIPT_PROP = 0.2
MAX_REF_TRANSCRIPT_PROP = 0.5
# width of the line that denote the gene or read
GENE_LINEWIDTH = 0.1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fusion-bed", "-f", type=Path, required=True, help=("The breakpoints groups file with the query fusion")
    )
    parser.add_argument(
        "--fusion-id", "-fid", type=str, required=True, help="Specify a fusion, such as BP600, to select from bed"
    )
    parser.add_argument("--bam", "-b", type=Path, help="Mapped bam file")
    parser.add_argument("--annotation", "-a", type=Path, help="Annotation gtf file")
    parser.add_argument("--pickle", "-p", type=Path, default=None, help="Use pickle file as input.")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Output png filename")

    args = parser.parse_args()

    return args


class BreakPoint(NamedTuple):
    """Defines the location of a break"""

    chrom: str
    pos: int
    is_reverse: bool
    gene_name: str
    gene_id: str


class BreakPointPair(NamedTuple):
    """Defines the two break points for a fusion event"""

    breakpoint1: BreakPoint
    breakpoint2: BreakPoint


class TranscriptFeature(str, Enum):
    EXON = "exon"
    CDS = "CDS"


class Interval(NamedTuple):
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


class TranscriptBlock(NamedTuple):
    coordinates: Interval
    feature: TranscriptFeature

    @property
    def start(self) -> int:
        return self.coordinates.start

    @property
    def end(self) -> int:
        return self.coordinates.end

    @property
    def length(self) -> int:
        return self.coordinates.length


class Transcript(NamedTuple):
    chrom: str
    blocks: List[TranscriptBlock]
    is_reverse: bool
    transcript_id: str
    gene_id: str

    @property
    def start(self) -> int:
        return min(block.start for block in self.blocks)

    @property
    def end(self) -> int:
        return max(block.end for block in self.blocks)


class SegmentBlocks(NamedTuple):
    """
    A simple storage to keep information from an alignment for plotting
    """

    chrom: str
    is_reverse: bool
    # a list of start and end positions of aligned gapless blocks.
    blocks: List[Interval]

    @property
    def start(self) -> int:
        return min(block.start for block in self.blocks)

    @property
    def end(self) -> int:
        return max(block.end for block in self.blocks)


@dataclass
class FusionAlignment:
    """A simple class that stores the "left" and "right" aligned segments."""

    read_name: str
    aligned_segment_5p: Optional[SegmentBlocks]
    aligned_segment_3p: Optional[SegmentBlocks]

    def __init__(self, read_name: str):
        self.read_name = read_name
        self.aligned_segment_5p = None
        self.aligned_segment_3p = None

    def add_alignment_blocks(self, chrom: str, is_reverse: bool, blocks: List[Interval], is_5p: bool):
        if is_5p:
            assert self.aligned_segment_5p is None
            self.aligned_segment_5p = SegmentBlocks(chrom, is_reverse, blocks)
        else:
            assert self.aligned_segment_3p is None
            self.aligned_segment_3p = SegmentBlocks(chrom, is_reverse, blocks)

    def get_5p_coordinate(self) -> int:
        """Returns the coordinate of the 5' extreme of the read"""
        if self.aligned_segment_5p.is_reverse:
            return self.aligned_segment_5p.end
        else:
            return self.aligned_segment_5p.start

    def get_3p_coordinate(self) -> int:
        """Returns the coordinate of the 3' extreme of the read"""
        if self.aligned_segment_3p.is_reverse:
            return self.aligned_segment_3p.start
        else:
            return self.aligned_segment_3p.end


def read_fusion_bedpe(
    bedpe_path: Path,
    fusion_id: str,
) -> Tuple[List[str], BreakPointPair]:
    """
    Parse the read names and breakpoints information for the specified fusion

    Args:
        bedpe_path: the path to the input file to parse
        fusion_id: the fusion ID to parse

    Returns:
        read_names: the list of read names provided in the file
        breakpoint_pair: the two breakpoints that form the fusion
    """
    breakpoint_pair = None
    with bedpe_path.open("rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            # make sure we have the correct number of fields
            line_fields = line.rstrip("\r\n").split("\t")
            if len(line_fields) != 12:
                raise ValueError(
                    f"Could not parse line {line.strip()} from bed PE file {bedpe_path}. Expected 11 fields but got "
                    f"{len(line_fields)} fields."
                )
            chr1, start1, _, chr2, start2, _, bp_id, _, strand1, strand2, info, extra = line_fields

            if bp_id != fusion_id:
                continue

            if breakpoint_pair is not None:
                raise ValueError(f"Found more than one entry for fusion ID {fusion_id}")

            # parse the read names from the extra field
            read_names = set()
            for field in extra.split(";"):
                if field.startswith("RN="):
                    read_names = set(field[3:].split(","))

            # parse the gene information from the info field
            gene_info = {}
            for field in info.split(";"):
                info_key, info_value = field.split("=")
                if info_key in ("GI", "GN"):
                    gene_info[info_key] = info_value.split(",")
                    if len(gene_info[info_key]) != 2:
                        raise ValueError(f"Expected two entries for info {info_key} on line {line.strip()}")
            for info_key in ("GI", "GN"):
                if info_key not in gene_info:
                    raise ValueError(f"No info tag {info_key} on line {line.strip()}")

            # create breakpoints
            breakpoint1 = BreakPoint(chr1, int(start1), strand1.strip() == "-", gene_info["GN"][0], gene_info["GI"][0])
            breakpoint2 = BreakPoint(chr2, int(start2), strand2.strip() == "-", gene_info["GN"][1], gene_info["GI"][1])
            breakpoint_pair = BreakPointPair(breakpoint1, breakpoint2)

    if breakpoint_pair is None:
        raise ValueError(f"Could not find entry for fusion ID {fusion_id}")

    return read_names, breakpoint_pair


def parse_transcripts_from_gtf(gtf_path: Path, gene_ids: List[str]) -> List[Transcript]:
    """Read the transcripts from the GTF file for the provided gene_ids

    Args:
        gtf_path: Input GTF file to parse
        gene_ids: The gene_ids for which to parse transcripts

    Returns:
        transcripts: the parsed transcripts
    """
    transcripts_by_id: Dict[str, Transcript] = {}
    with gtf_path.open("rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            line_fields = line.strip().split("\t")
            if len(line_fields) != 9:
                raise ValueError(
                    f"Could not parse line {line.strip()} from GTF file {gtf_path}. Expected 9 fields but got "
                    f"{len(line_fields)} fields."
                )
            chrom, source, feature, start, end, score, strand, frame, attributes = line_fields

            if feature in ["exon", "CDS"]:
                start = int(start)
                end = int(end)
                transcript_id = None
                gene_id = None
                # the last attribute ends with a `;` so we need to ignore the last field
                for attribute in attributes.split(";")[:-1]:
                    try:
                        key, value = attribute.strip().split(maxsplit=1)
                    except Exception as exc:
                        raise ValueError(f"Attribute {attribute} has improper format. Line {line.strip()}") from exc
                    if key == "transcript_id":
                        transcript_id = value.strip('"')
                    elif key == "gene_id":
                        gene_id = value.strip('"')
                if gene_id is None or transcript_id is None:
                    raise ValueError(f"Missing gene_id or transcript_id from GTF on line {line.strip()}")
                if gene_id not in gene_ids:
                    continue

                is_reverse = strand == "-"
                transcript_block = TranscriptBlock(coordinates=Interval(start, end), feature=TranscriptFeature(feature))

                if transcript_id in transcripts_by_id:
                    # check that we have the same chromosome, strand and gene_id
                    assert chrom == transcripts_by_id[transcript_id].chrom
                    assert gene_id == transcripts_by_id[transcript_id].gene_id
                    assert is_reverse == transcripts_by_id[transcript_id].is_reverse
                    # and add the block
                    transcripts_by_id[transcript_id].blocks.append(transcript_block)
                else:
                    # store the new transcript
                    transcripts_by_id[transcript_id] = Transcript(
                        chrom=chrom,
                        blocks=[transcript_block],
                        is_reverse=is_reverse,
                        transcript_id=transcript_id,
                        gene_id=gene_id,
                    )

    return list(transcripts_by_id.values())


def read_fusion_alignments(bam_path: Path, read_names: List[str]) -> List[FusionAlignment]:
    """Get alignments from bam for the reads supporting the fusion. We ensure to keep all the blocks of the alignment.

    Args:
        bam_path: input bam
        read_names: fusion supporting query names

    Returns:
        alignments: fusion alignments
    """
    pysam.set_verbosity(0)
    open_mode = "r" if bam_path.suffix == "sam" else "rb"

    read_names = set(read_names)
    alignments_by_name: Dict[str, FusionAlignment] = {}
    with pysam.AlignmentFile(str(bam_path), open_mode, check_sq=False, require_index=False) as bam_file:
        for alignment in bam_file.fetch(until_eof=True):
            if alignment.query_name in read_names:
                # TODO: not sure why we're excluding reads with more that 2 suppl alignments
                if len(alignment.get_tag("SA").split(";")) > 2:
                    continue

                # check if this alignment is the 5p end of the molecule
                # TODO: we should get rid of this magic "50 bps from the start of the read"
                #    recommendation would be to store the blocks for all aligned segment for our read, with the
                #    start coordinate of the read.
                is_5p = (not alignment.is_reverse and alignment.query_alignment_start < 50) or (
                    not alignment.is_reverse and alignment.query_length - alignment.query_alignment_end < 50
                )
                blocks = [Interval(*block) for block in alignment.get_blocks()]

                if alignment.query_name not in alignments_by_name:
                    alignments_by_name[alignment.query_name] = FusionAlignment(alignment.query_name)
                try:
                    alignments_by_name[alignment.query_name].add_alignment_blocks(
                        chrom=alignment.reference_name, is_reverse=alignment.is_reverse, blocks=blocks, is_5p=is_5p
                    )
                except Exception as exc:
                    raise ValueError(
                        f"Error adding blocks for {alignment.query_name} at "
                        f"{alignment.reference_name}:{alignment.query_alignment_start}"
                    ) from exc

    return list(alignments_by_name.values())


def filter_alignments(
    fusion_alignments: List[FusionAlignment], is_reverse_5p: bool, is_reverse_3p: bool
) -> List[FusionAlignment]:
    """
    Filter alignments and keep those that have both sides defined, and the correct strand on each side.
    Log warnings for any omitted alignment.
    """

    def keep_alignment(alignment: FusionAlignment):
        if alignment.aligned_segment_5p is None or alignment.aligned_segment_3p is None:
            logger.warning(f"Alignment for {alignment.read_name} is missing one side. Omitting from plot.")
            return False
        if alignment.aligned_segment_5p.is_reverse != is_reverse_5p:
            logger.warning(f"5p alignment for {alignment.read_name} is on the wrong strand. Omitting from plot.")
            return False
        if alignment.aligned_segment_3p.is_reverse != is_reverse_3p:
            logger.warning(f"3p alignment for {alignment.read_name} is on the wrong strand. Omitting from plot.")
            return False
        return True

    return list(filter(keep_alignment, fusion_alignments))


def filter_ref_transcripts(
    ref_transcripts: List[Transcript],
    gene_id: str,
    is_reverse: bool,
    overlap: Interval,
) -> List[Transcript]:
    """
    Filter transcripts for those on the provided gene, the correct strand and some overlap with our region of interest.
    Log warning for wrong strand.
    """

    def keep_transcript(transcript: Transcript):
        if transcript.gene_id != gene_id:
            return False
        if (
            max(block.end for block in transcript.blocks) < overlap.start
            or min(block.start for block in transcript.blocks) > overlap.end
        ):
            return False
        if transcript.is_reverse != is_reverse:
            logger.warning(f"Transcript {transcript.transcript_id} is on the wrong strand. Omitting from plot.")
            return False
        return True

    return list(filter(keep_transcript, ref_transcripts))


def plot_fusion(
    ref_transcripts: List[Transcript],
    fusion_alignments: List[FusionAlignment],
    breakpoints: BreakPointPair,
) -> plt.Figure:
    """Create a fusion figure that shows the genes annotation and reads that support the fusion event.

    Args:
        ref_transcripts: transcripts from the reference annotation
        fusion_alignments: the alignments for the reads supporting the fusion
        breakpoints: the breakpoints for the called fusion

    Returns:
        figure: a matplotlib figure depicting the fusion event
    """
    breakpoint1, breakpoint2 = breakpoints

    # filter the alignments
    alignments: List[FusionAlignment] = filter_alignments(
        fusion_alignments, breakpoint1.is_reverse, breakpoint2.is_reverse
    )
    if len(alignments) == 0:
        raise ValueError("Could not find any alignments for provided fusion.")

    # get the extreme coordinates across the alignments
    gene1_start = (
        max(fa.get_5p_coordinate() for fa in alignments)
        if breakpoint1.is_reverse
        else min(fa.get_5p_coordinate() for fa in alignments)
    )
    gene2_end = (
        min(fa.get_3p_coordinate() for fa in alignments)
        if breakpoint2.is_reverse
        else max(fa.get_3p_coordinate() for fa in alignments)
    )

    # filter the transcripts
    gene1_transcripts: List[Transcript] = filter_ref_transcripts(
        ref_transcripts, breakpoint1.gene_id, breakpoint1.is_reverse, Interval(*sorted((gene1_start, breakpoint1.pos)))
    )
    gene2_transcripts: List[Transcript] = filter_ref_transcripts(
        ref_transcripts, breakpoint2.gene_id, breakpoint2.is_reverse, Interval(*sorted((gene2_end, breakpoint2.pos)))
    )

    # compute the span we'll need for each side of the fusion
    gene1_start = (
        max(fa.get_5p_coordinate() for fa in alignments)
        if breakpoint1.is_reverse
        else min(fa.get_5p_coordinate() for fa in alignments)
    )
    gene2_end = (
        min(fa.get_3p_coordinate() for fa in alignments)
        if breakpoint2.is_reverse
        else max(fa.get_3p_coordinate() for fa in alignments)
    )
    gene1_span = abs(breakpoint1.pos - gene1_start)
    gene2_span = abs(breakpoint2.pos - gene2_end)
    gene1_xprop = gene1_span / (gene1_span + gene2_span)
    gene2_xprop = 1 - gene1_xprop

    # create the figures and panels according to the computed proportions
    total_plot_transcripts = len(gene1_transcripts) + len(gene2_transcripts)
    plot_prop = 1.0 - 2 * PLOT_MARGIN
    ref_transcript_prop = total_plot_transcripts / (total_plot_transcripts + len(alignments))
    ref_transcript_prop = min(MAX_REF_TRANSCRIPT_PROP, ref_transcript_prop)
    ref_transcript_prop = max(MIN_REF_TRANSCRIPT_PROP, ref_transcript_prop)
    reads_yprop = (1.0 - ref_transcript_prop) * plot_prop
    # assigning gene1 / gene2 props based on number of transcripts for each gene.
    gene1_yprop = len(gene1_transcripts) / total_plot_transcripts * ref_transcript_prop * plot_prop
    gene2_yprop = len(gene2_transcripts) / total_plot_transcripts * ref_transcript_prop * plot_prop

    figure = plt.figure(figsize=(8, 3))
    gene1_panel = plt.axes([PLOT_MARGIN, PLOT_MARGIN + reads_yprop + gene2_yprop, plot_prop, gene1_yprop])
    gene2_panel = plt.axes([PLOT_MARGIN, PLOT_MARGIN + reads_yprop, plot_prop, gene2_yprop])
    reads1_panel = plt.axes([PLOT_MARGIN, PLOT_MARGIN, gene1_xprop * plot_prop, reads_yprop])
    reads2_panel = plt.axes(
        [PLOT_MARGIN + (gene1_xprop * plot_prop), PLOT_MARGIN, (gene2_xprop * plot_prop), reads_yprop]
    )

    # plot the gene annotations
    for panel, transcripts, color in (
        (gene1_panel, gene1_transcripts, GENE1_COLOR),
        (gene2_panel, gene2_transcripts, GENE2_COLOR),
    ):
        for ypos, transcript in enumerate(transcripts, start=1):
            panel.plot([transcript.start, transcript.end], [ypos] * 2, lw=GENE_LINEWIDTH, c=color, zorder=10)
            for transcript_block in transcript.blocks:
                offset = 0.125 if transcript_block.feature == TranscriptFeature.EXON else 0.25
                panel.add_patch(
                    Rect(
                        (transcript_block.coordinates.start, ypos - offset),
                        transcript_block.coordinates.length,
                        offset * 2,
                        lw=0,
                        fc=color,
                        zorder=12,
                    )
                )

    # plot the reads
    for ypos, falignment in enumerate(
        sorted(alignments, key=lambda fa: fa.get_5p_coordinate(), reverse=breakpoint1.is_reverse), start=1
    ):
        # plot the 5' side of the fusion
        for block in falignment.aligned_segment_5p.blocks:
            reads1_panel.add_patch(Rect((block.start, ypos - 0.25), block.length, 0.5, lw=0, fc=GENE1_COLOR, zorder=12))
        reads1_panel.plot(
            [falignment.get_5p_coordinate(), breakpoint1.pos], [ypos] * 2, lw=GENE_LINEWIDTH, c=GENE1_COLOR, zorder=10
        )

        # plot the 3' side of the fusion
        for block in falignment.aligned_segment_3p.blocks:
            reads2_panel.add_patch(Rect((block.start, ypos - 0.25), block.length, 0.5, lw=0, fc=GENE2_COLOR, zorder=12))
        reads2_panel.plot(
            [falignment.get_3p_coordinate(), breakpoint2.pos], [ypos] * 2, lw=GENE_LINEWIDTH, c=GENE2_COLOR, zorder=10
        )

    gene1_panel.set_ylabel(f"{breakpoint1.gene_name}\n({'-' if breakpoint1.is_reverse else '+'})", fontsize=5)
    gene2_panel.set_ylabel(f"{breakpoint2.gene_name}\n({'-' if breakpoint2.is_reverse else '+'})", fontsize=5)
    reads1_panel.set_ylabel("Reads")

    # set the limits on each panel
    if breakpoint1.is_reverse:
        # this will invert the axis
        assert gene1_start > breakpoint1.pos
        gene1_panel.set_xlim(gene1_start, breakpoint1.pos - gene2_span)
    else:
        assert gene1_start < breakpoint1.pos
        gene1_panel.set_xlim(gene1_start, breakpoint1.pos + gene2_span)
    if breakpoint2.is_reverse:
        # this will invert the axis
        assert gene2_end < breakpoint2.pos
        gene2_panel.set_xlim(breakpoint2.pos + gene1_span, gene2_end)
    else:
        assert gene2_end > breakpoint2.pos
        gene2_panel.set_xlim(breakpoint2.pos - gene1_span, gene2_end)
    # Note: if either of these are on the - strand, that axis will get inverted.
    reads1_panel.set_xlim(gene1_start, breakpoint1.pos)
    reads2_panel.set_xlim(breakpoint2.pos, gene2_end)

    reads1_panel.set_ylim(0, len(alignments) + 1)
    reads2_panel.set_ylim(0, len(alignments) + 1)

    reads1_panel.spines["right"].set_visible(False)
    reads2_panel.spines["left"].set_visible(False)
    for panel in [gene1_panel, gene2_panel, reads1_panel, reads2_panel]:
        panel.tick_params(
            bottom=False,
            labelbottom=False,
            left=False,
            labelleft=False,
            right=False,
            labelright=False,
            top=False,
            labeltop=False,
        )
        figure.add_axes(panel)
    return figure


def main(args: argparse.Namespace) -> int:

    try:
        if args.pickle is None and (args.bam is None or args.annotation is None):
            raise ValueError("Please specify input. This can be a pickle file, or a bam and annotation")
        pickle_opath = Path(f"pbfusion_viz.{args.fusion_id}.pickle")
        read_names, breakpoint_pair = read_fusion_bedpe(args.fusion_bed, args.fusion_id)
        if args.pickle is not None:
            with args.pickle.open("rb") as pickle_stream:
                ref_transcripts = pickle.load(pickle_stream)
                alignments = pickle.load(pickle_stream)
        else:
            ref_transcripts = parse_transcripts_from_gtf(
                args.annotation, [breakpoint.gene_id for breakpoint in breakpoint_pair]
            )
            alignments = read_fusion_alignments(args.bam, read_names)
            with pickle_opath.open("wb") as pickle_stream:
                pickle.dump(ref_transcripts, pickle_stream)
                pickle.dump(alignments, pickle_stream)

        figure = plot_fusion(ref_transcripts, alignments, breakpoint_pair)
        output_png = args.output if args.output.suffix == ".png" else args.output.with_suffix(".png")
        figure.savefig(output_png, dpi=1000)
    except Exception as exc:
        logger.exception(f"FATAL error encountered {exc}")

    return 0


if __name__ == "__main__":
    args = parse_args()
    sys.exit(main(args))
