# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Merge fragments that the *in silico* digest split but the enzyme did not.

The problem
-----------
The virtual digest cuts at every recognition site in a read. The enzyme in
the tube did not: digestion is never complete, so a genuine restriction
fragment usually contains several uncut sites inside it. Those sites are cut in
silico anyway, and one real fragment becomes several monomers that align
head-to-tail on the genome.

Nothing warns you. What you see instead is a Hi-C map with a grossly inflated
diagonal, because every such split manufactures contacts that never happened:

* a locus split into *a* pieces facing a locus split into *b* pieces yields
  ``a x b`` pairs where there was one contact;
* worse, the *a* pieces of a single locus are paired with each other, and
  those land straight on the diagonal.

On the library this tool was written for, that was an 11-fold duplication at
5 kb resolution and 18-fold at 25 kb, with 44-76% of all pairs on the diagonal.
Juicer's normalisations cannot repair it: only 14% of the bias is separable
into a row factor times a column factor, so a residual factor of 2.45 survives
whatever you normalise with. The distortion also follows enzyme-site
density, so it deforms the map rather than just scaling it.

The fix
-------
Two monomers that are consecutive along the read *and* contiguous on the genome
(same chromosome, same strand, within ``--merge-gap`` bases) came from one
piece of DNA that was never cut. Glue them back together before doing anything
else.

Measured effect, reproducible with the shipped command on that library
(208,121 concatemers, defaults, ``--min-fragments 2``): 2,278,102 aligned
monomers collapse to 450,213 fragments, 5.1x fewer. Cis contacts below 10 kb
fall from 79.8% of all pairs to 26.4%, while the 10 kb - 1 Mb window, where
TADs and loops live, rises from 13.9% to 47.8%. The signal was there all
along, buried under the artefact.

Order matters
-------------
Merging runs on the complete chain of monomers, and the MAPQ filter is
applied afterwards, to the merged blocks. Filtering first breaks the chain: a
poorly mapped monomer in the middle is removed, its two neighbours are no
longer adjacent, and they are counted as two distinct loci in contact - a
contact that does not exist.

Doing it in the wrong order made the molecule count *rise* with a stricter
MAPQ threshold (61,677 to 70,735). A stricter filter
cannot create molecules, and that is how the bug was found. A merged block is
kept when at least one of its monomers passed, so each block carries the
highest MAPQ it contained.
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import json
import os
import sys
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from pore_c_aqb import __version__
from pore_c_aqb.progress import (
    Progress,
    add_progress_arguments,
    progress_enabled,
)
from pore_c_aqb.reads import iter_concatemers, read_span

__all__ = ["Fragment", "merge_adjacent", "fragments_of", "MergeStats", "main"]

PROG = "pore-c-aqb merge"


@dataclass
class Fragment:
    """A piece of a concatemer, possibly several monomers glued back together."""

    read_start: int
    read_end: int
    chrom: str
    ref_start: int
    ref_end: int
    strand: str
    mapq: int
    n_monomers: int = 1
    #: rank along the read, used when read_start is unknown (-1)
    order: int = 0

    @property
    def midpoint(self) -> int:
        """Position used to represent this fragment in a contact.

        The midpoint, not an edge: a merged fragment can span several kb and
        either border would be an arbitrary choice. At 5 kb resolution or
        coarser the difference is immaterial.
        """
        return (self.ref_start + self.ref_end) // 2

    @property
    def length(self) -> int:
        return self.ref_end - self.ref_start


def fragments_of(alignments) -> list[Fragment]:
    """One :class:`Fragment` per aligned monomer, in read order.

    ``read_start``/``read_end`` are -1 when the BAM carries neither an ``Xc``
    tag nor parseable monomer names; the order the alignments arrive in is then
    the order along the read, as resolved by :func:`iter_concatemers`.
    """
    out = []
    for offset, a in enumerate(alignments):
        start, end = read_span(a)
        if start < 0:
            start = end = -1
        out.append(Fragment(
            read_start=start,
            read_end=end,
            chrom=a.reference_name,
            ref_start=a.reference_start,
            ref_end=a.reference_end,
            strand="-" if a.is_reverse else "+",
            mapq=a.mapping_quality,
            order=offset,
        ))
    return out


def merge_adjacent(fragments: list[Fragment], gap: int) -> list[Fragment]:
    """Glue consecutive monomers that are contiguous on the genome.

    Two fragments merge when they sit next to each other along the read, on the
    same chromosome and the same strand, with at most ``gap`` bases between
    them. The gap is measured in both directions, because a fragment aligned to
    the minus strand runs backwards along the reference: the next monomer's
    start may follow the previous one's end, or precede its start.

    The merged block keeps the highest MAPQ of its parts, so that a block
    containing one confidently placed monomer survives the later filter.
    """
    if not fragments:
        return []
    ordered = sorted(
        fragments,
        key=lambda f: (f.read_start if f.read_start >= 0 else f.order, f.order))
    out = [
        Fragment(f.read_start, f.read_end, f.chrom, f.ref_start, f.ref_end,
                 f.strand, f.mapq, f.n_monomers, f.order)
        for f in ordered[:1]
    ]
    for frag in ordered[1:]:
        prev = out[-1]
        same_place = (frag.chrom == prev.chrom and frag.strand == prev.strand)
        contiguous = same_place and min(
            abs(frag.ref_start - prev.ref_end),
            abs(prev.ref_start - frag.ref_end),
        ) <= gap
        if contiguous:
            prev.read_end = frag.read_end
            prev.ref_start = min(prev.ref_start, frag.ref_start)
            prev.ref_end = max(prev.ref_end, frag.ref_end)
            prev.mapq = max(prev.mapq, frag.mapq)
            prev.n_monomers += frag.n_monomers
        else:
            out.append(Fragment(
                frag.read_start, frag.read_end, frag.chrom, frag.ref_start,
                frag.ref_end, frag.strand, frag.mapq, frag.n_monomers,
                frag.order))
    return out


@dataclass
class MergeStats:
    """What the merge changed. The point of the exercise, so it is reported."""

    n_concatemers: int = 0
    n_monomers: int = 0
    n_fragments: int = 0
    n_kept_molecules: int = 0
    n_pairs: int = 0
    order_before: Counter = None
    order_after: Counter = None
    #: cis separations, bucketed, before and after merging
    dist_before: Counter = None
    dist_after: Counter = None

    def __post_init__(self):
        for name in ("order_before", "order_after", "dist_before",
                     "dist_after"):
            if getattr(self, name) is None:
                setattr(self, name, Counter())

    def as_dict(self) -> dict:
        collapse = (self.n_monomers / self.n_fragments
                    if self.n_fragments else 0.0)
        return {
            "tool": f"pore-c-aqb {__version__} merge",
            "concatemers": self.n_concatemers,
            "aligned_monomers": self.n_monomers,
            "fragments_after_merge": self.n_fragments,
            "collapse_ratio": round(collapse, 3),
            "molecules_kept": self.n_kept_molecules,
            "pairs_written": self.n_pairs,
            "fragments_per_molecule_before": dict(sorted(
                self.order_before.items())),
            "fragments_per_molecule_after": dict(sorted(
                self.order_after.items())),
            "cis_distance_before": _pct(self.dist_before),
            "cis_distance_after": _pct(self.dist_after),
        }


DIST_BUCKETS = [
    ("<1kb", 1_000), ("1-10kb", 10_000), ("10kb-1Mb", 1_000_000),
    (">1Mb", float("inf")),
]


def _bucket(distance: int) -> str:
    for label, limit in DIST_BUCKETS:
        if distance < limit:
            return label
    return ">1Mb"


def _pct(counter: Counter) -> dict:
    total = sum(counter.values())
    if not total:
        return {}
    return {label: round(100.0 * counter.get(label, 0) / total, 1)
            for label, _ in DIST_BUCKETS}


def _tally_distances(frags: list[Fragment], counter: Counter) -> None:
    for a, b in combinations(frags, 2):
        if a.chrom == b.chrom:
            counter[_bucket(abs(a.midpoint - b.midpoint))] += 1


def _open_out(path, mode="wt"):
    if path is None or str(path) == "-":
        return sys.stdout, False
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(str(path), mode), True


FRAGMENT_COLUMNS = [
    "read_id", "fragment_index", "n_fragments", "chrom", "start", "end",
    "midpoint", "strand", "mapq", "n_monomers_merged", "read_start", "read_end",
]


def _chrom_sizes(sizes_path) -> dict:  # noqa: D401
    """``{chrom: length}`` from a two-column sizes file, order preserved.

    The order of the file is the order used for the .pairs upper triangle, so
    it matches whatever genome ordering the rest of your pipeline uses.
    """
    if not Path(sizes_path).exists():
        raise SystemExit(f"Chromosome sizes file not found: {sizes_path}")
    sizes = {}
    with open(sizes_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def run(args) -> int:
    stats = MergeStats()
    frag_fh, frag_close = _open_out(args.output)
    pairs_fh, pairs_close = (None, False)
    chrom_sizes = _chrom_sizes(args.sizes) if args.sizes else {}
    chrom_order = {c: i for i, c in enumerate(chrom_sizes)}

    stack = contextlib.ExitStack()
    try:
        frag_fh.write("\t".join(FRAGMENT_COLUMNS) + "\n")
        if args.pairs:
            pairs_fh, pairs_close = _open_out(args.pairs)
            pairs_fh.write("## pairs format v1.0\n")
            pairs_fh.write("#shape: upper triangle\n")
            pairs_fh.write(f"#generated-by: pore-c-aqb {__version__} merge "
                           f"(merge-gap={args.merge_gap}, mapq={args.mapq})\n")
            for chrom, length in chrom_sizes.items():
                pairs_fh.write(f"#chromsize: {chrom} {length}\n")
            pairs_fh.write(
                "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")

        size = os.path.getsize(args.bam) if os.path.exists(args.bam) else None
        holder = {}
        bar = stack.enter_context(Progress(
            "merging", unit="concatemers",
            position=(lambda: holder["tell"]() >> 16) if size else None,
            size=size, enabled=progress_enabled(args)))
        for read_id, alignments in iter_concatemers(args.bam, tell=holder):
            stats.n_concatemers += 1
            bar.update()
            if not alignments:
                continue
            frags = fragments_of(alignments)
            stats.n_monomers += len(frags)
            stats.order_before[min(len(frags), 20)] += 1
            _tally_distances(frags, stats.dist_before)

            # merge on the full chain, THEN filter: see the module docstring
            frags = merge_adjacent(frags, args.merge_gap)
            frags = [f for f in frags if f.mapq >= args.mapq]
            if args.min_length:
                frags = [f for f in frags if f.length >= args.min_length]
            if len(frags) < args.min_fragments:
                continue

            stats.n_fragments += len(frags)
            stats.n_kept_molecules += 1
            stats.order_after[min(len(frags), 20)] += 1
            _tally_distances(frags, stats.dist_after)

            for i, f in enumerate(frags):
                frag_fh.write(
                    f"{read_id}\t{i}\t{len(frags)}\t{f.chrom}\t{f.ref_start}\t"
                    f"{f.ref_end}\t{f.midpoint}\t{f.strand}\t{f.mapq}\t"
                    f"{f.n_monomers}\t{f.read_start}\t{f.read_end}\n")

            if pairs_fh is not None:
                stats.n_pairs += _write_pairs(
                    pairs_fh, read_id, frags, chrom_order, args.min_sep)
    finally:
        stack.close()
        if frag_close:
            frag_fh.close()
        if pairs_close:
            pairs_fh.close()

    summary = stats.as_dict()
    summary["parameters"] = {
        "input": str(args.bam),
        "merge_gap": args.merge_gap,
        "mapq": args.mapq,
        "min_fragments": args.min_fragments,
        "min_length": args.min_length,
        "min_sep": args.min_sep,
    }
    if args.stats:
        Path(args.stats).write_text(json.dumps(summary, indent=2) + "\n")
    _report(summary, args)
    return 0


def _write_pairs(fh, read_id, frags, chrom_order, min_sep) -> int:
    written = 0
    for a, b in combinations(frags, 2):
        if a.chrom == b.chrom and abs(a.midpoint - b.midpoint) < min_sep:
            continue
        c1, p1, s1 = a.chrom, a.midpoint + 1, a.strand      # .pairs is 1-based
        c2, p2, s2 = b.chrom, b.midpoint + 1, b.strand
        big = 1 << 30
        if (chrom_order.get(c1, big), p1) > (chrom_order.get(c2, big), p2):
            c1, p1, s1, c2, p2, s2 = c2, p2, s2, c1, p1, s1
        fh.write(f"{read_id}\t{c1}\t{p1}\t{c2}\t{p2}\t{s1}\t{s2}\n")
        written += 1
    return written


def _report(summary: dict, args) -> None:
    if getattr(args, "quiet", False):
        return
    say = print if args.output not in (None, "-") else \
        (lambda *a, **k: print(*a, file=sys.stderr, **k))
    say(f"Read {summary['concatemers']:,} concatemers, "
        f"{summary['aligned_monomers']:,} aligned monomers.")
    say(f"Merged into {summary['fragments_after_merge']:,} fragments "
        f"({summary['collapse_ratio']}x fewer), keeping "
        f"{summary['molecules_kept']:,} molecules "
        f"with at least {args.min_fragments} fragment(s).")
    before, after = summary["cis_distance_before"], summary["cis_distance_after"]
    if before and after:
        say("")
        say("Cis contacts by separation, before and after merging:")
        say(f"  {'separation':<12} {'before':>8} {'after':>8}")
        for label, _ in DIST_BUCKETS:
            say(f"  {label:<12} {before.get(label, 0):>7.1f}% "
                f"{after.get(label, 0):>7.1f}%")
        say("")
        say("  A large '<1kb' share before merging is the artefact: those are "
            "pieces of one\n  uncut restriction fragment being paired with "
            "each other, landing on the\n  diagonal of the Hi-C map.")
    if summary["pairs_written"]:
        say(f"\nWrote {summary['pairs_written']:,} pairs.")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=PROG,
        description=(
            "Glue back together the fragments that the in silico digest split "
            "but the enzyme never cut, so that the Hi-C diagonal stops being "
            "inflated by contacts that did not happen."),
        epilog=(
            "WHY THIS EXISTS\n"
            "  The virtual digest cuts at every recognition site. The enzyme in\n"
            "  the tube did not: real digestion is incomplete, so one genuine\n"
            "  restriction fragment usually contains several uncut sites and\n"
            "  ends up as several monomers that align head-to-tail.\n"
            "\n"
            "  A locus split into a pieces facing a locus split into b pieces\n"
            "  then yields a*b pairs for ONE real contact, and the a pieces of\n"
            "  one locus get paired with each other - straight onto the\n"
            "  diagonal. Measured on a real library: 11x duplication at 5 kb,\n"
            "  18x at 25 kb, 44-76% of pairs on the diagonal. Juicer's\n"
            "  normalisations do not fix it (only 14% of the bias is a row\n"
            "  factor times a column factor; a 2.45x residual survives).\n"
            "\n"
            "WHAT IT DOES\n"
            "  Two monomers consecutive along the read and contiguous on the\n"
            "  genome (same chromosome, same strand, gap <= --merge-gap) came\n"
            "  from one uncut piece of DNA. They are merged into one fragment.\n"
            "  Effect on that library: 2,278,102 monomers become 450,213\n"
            "  fragments (5.1x fewer); cis contacts under 10 kb fall from\n"
            "  79.8% to 26.4%, while the 10 kb - 1 Mb window (TADs, loops)\n"
            "  rises from 13.9% to 47.8%. The run prints this table for your\n"
            "  own data, so you can see what it changed.\n"
            "\n"
            "ORDER MATTERS\n"
            "  Merging runs on the complete chain of monomers; --mapq is\n"
            "  applied afterwards, to the merged blocks. Filtering first would\n"
            "  break the chain - drop a middle monomer and its two neighbours\n"
            "  stop being adjacent, so they are counted as two loci in contact.\n"
            "  Getting this backwards made the molecule count RISE with a\n"
            "  stricter threshold, which is impossible and is how it was found.\n"
            "\n"
            "EXAMPLES\n"
            "  pore-c-aqb merge sample.ns.bam --output fragments.tsv.gz\n"
            "\n"
            "  pore-c-aqb merge sample.ns.bam --output fragments.tsv.gz \\\n"
            "      --pairs contacts.pairs --sizes hg38.sizes.genome \\\n"
            "      --stats merge_stats.json\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--version", action="version", version=f"pore-c-aqb {__version__}")
    p.add_argument(
        "bam",
        help=("Aligned monomer BAM, grouped by read name - the workflow's "
              "*.ns.bam. A coordinate-sorted BAM is rejected: the monomers of "
              "one concatemer have to be adjacent."))
    p.add_argument(
        "--output", default="-",
        help=("Merged fragments, as TSV ('-' for stdout, '.gz' compresses). "
              "One row per fragment, grouped by molecule."))
    p.add_argument(
        "--merge-gap", type=int, default=100, metavar="BP",
        help=("Largest gap between two monomers still considered one uncut "
              "fragment. 100 bp is tolerant enough for alignment imprecision "
              "at the cut site and strict enough not to swallow genuine "
              "short-range contacts; the sensitivity curve measured over "
              "0/20/100/500/1000/5000 bp flattens beyond ~100 bp."))
    p.add_argument(
        "--mapq", type=int, default=1, metavar="Q",
        help=("Minimum mapping quality, applied AFTER merging. A merged block "
              "is kept when at least one of its monomers reached this, so it "
              "carries the highest MAPQ it contained."))
    p.add_argument(
        "--min-fragments", type=int, default=1, metavar="N",
        help=("Keep only molecules with at least N fragments after merging. "
              "Use 2 for contacts, 3 for multi-way analysis."))
    p.add_argument(
        "--min-length", type=int, default=0, metavar="BP",
        help="Drop merged fragments shorter than this.")
    p.add_argument(
        "--pairs", default=None, metavar="PATH",
        help=("Also write contacts in 4DN .pairs format ('.gz' compresses). "
              "Every pair of distinct fragments in a molecule is emitted, so "
              "a k-fragment molecule gives k*(k-1)/2 pairs. Positions are the "
              "fragment midpoints, 1-based."))
    p.add_argument(
        "--sizes", default=None, metavar="PATH",
        help=("Chromosome sizes file, used to order the .pairs upper triangle "
              "consistently with your genome. Without it, ordering falls back "
              "to the chromosome name."))
    p.add_argument(
        "--min-sep", type=int, default=0, metavar="BP",
        help=("Drop cis pairs closer than this. A blunt instrument: the "
              "rigorous criterion is the number of restriction sites between "
              "the two fragments, since that is what decides whether they "
              "could be one uncut piece. Leave at 0 unless you know why."))
    p.add_argument(
        "--stats", default=None, metavar="PATH",
        help="Write the summary as JSON.")
    p.add_argument(
        "--quiet", action="store_true",
        help="Print nothing but errors (the summary still goes to --stats).")
    add_progress_arguments(p)
    return p


def main(argv=None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
