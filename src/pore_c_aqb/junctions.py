#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Which enzymes actually cut? Ask the aligned data.

The digest report counts recognition sites present in the reads. That number
cannot tell a working enzyme from one that never left the freezer, because
every motif occurs in genomic DNA by chance. This script performs the test that
can, and it needs alignments.

The idea
--------
Inside a concatemer, the boundary between two consecutive monomers that land in
different places in the genome is a **ligation junction**: two ends that were
cut by an enzyme and joined back together. If enzyme E cut, its site sits at
those boundaries far more often than at random positions. If E never cut, its
site is at exactly background frequency there.

So: collect junction boundaries, ask how often each enzyme's site is at one,
and divide by the same quantity measured at random positions on the same
chromosomes. That ratio is the answer. Around 1 means the enzyme did nothing.

The outer ends of a concatemer are deliberately ignored — they carry adapters
and trimming artefacts, not clean cut sites.

Usage
-----
    pore-c-aqb-junctions monomers.aligned.ns.bam ref.fa \\
        --enzymes DpnII,NlaIII,HindIII

The BAM must be name-sorted (as produced by the workflow) so that monomers of
one concatemer are adjacent, and must carry the ``MI`` concatemer tag. The
reference FASTA needs its ``.fai`` index.
"""
from __future__ import annotations

import argparse
import random
import sys
from collections import defaultdict

import pysam

from pore_c_aqb.enzymes import EnzymeSpecError, resolve_enzymes
from pore_c_aqb.sites import site_regex


def read_order_key(aln) -> int:
    """Position of this monomer along its concatemer.

    ``Xc`` is the (start, end, index, count) tag written by the digest. Falling
    back to the name suffix keeps the script usable on BAMs where Xc was
    stripped, since monomer names are ``<read>:<start>:<end>``.
    """
    try:
        return int(aln.get_tag("Xc")[0])
    except (KeyError, TypeError, IndexError):
        parts = aln.query_name.rsplit(":", 2)
        return int(parts[1]) if len(parts) == 3 else 0


def junction_boundaries(bam_path: str, mapq: int, min_jump: int):
    """Plus-strand coordinates of internal ligation junctions.

    Yields ``(chrom, position)`` for both sides of each junction where the two
    monomers are on different chromosomes or more than ``min_jump`` apart.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    groups: dict[str, list] = defaultdict(list)
    current = None

    def flush(alns):
        alns = [a for a in alns
                if not a.is_unmapped and not a.is_secondary
                and not a.is_supplementary and a.mapping_quality >= mapq]
        alns.sort(key=read_order_key)
        for left, right in zip(alns, alns[1:]):
            far = (left.reference_name != right.reference_name
                   or abs(right.reference_start - left.reference_start)
                   > min_jump)
            if not far:
                continue
            # the read-3' end of the left monomer and the read-5' end of the
            # right one; on a reverse alignment those swap in the genome
            yield (left.reference_name,
                   left.reference_start if left.is_reverse
                   else left.reference_end)
            yield (right.reference_name,
                   right.reference_end if right.is_reverse
                   else right.reference_start)

    for aln in bam:
        mi = aln.get_tag("MI") if aln.has_tag("MI") else \
            aln.query_name.rsplit(":", 2)[0]
        if current is not None and mi != current:
            yield from flush(groups.pop(current, []))
        current = mi
        groups[mi].append(aln)
    if current is not None:
        yield from flush(groups.pop(current, []))
    bam.close()


def motif_at(ref, chrom: str, pos: int, enzyme, tol: int) -> bool:
    """Is ``enzyme``'s site placed so that its cut falls within ``tol`` of pos?"""
    site, fst5, n = enzyme.site, enzyme.fst5, len(enzyme.site)
    lo = max(0, pos - fst5 - tol)
    hi = pos - fst5 + tol + n
    try:
        window = ref.fetch(chrom, lo, hi).upper()
    except (ValueError, KeyError):
        return False
    return site_regex(site).search(window) is not None


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__.split("Usage")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam", help="Name-sorted aligned monomer BAM.")
    p.add_argument("reference", help="Indexed reference FASTA (.fai needed).")
    p.add_argument("--enzymes", required=True,
                   help="Comma-separated names, e.g. 'DpnII,NlaIII,HindIII'.")
    p.add_argument("--mapq", type=int, default=20,
                   help="Minimum monomer MAPQ.")
    p.add_argument("--min-jump", type=int, default=1000,
                   help="Minimum genomic distance for a boundary to count as "
                        "a ligation junction rather than an uncut site.")
    p.add_argument("--tol", type=int, default=2,
                   help="Slack in bases between the motif's cut and the "
                        "alignment boundary; covers the overhang.")
    p.add_argument("--max-junctions", type=int, default=200000,
                   help="Stop after this many boundaries (0 = all).")
    p.add_argument("--seed", type=int, default=0,
                   help="Seed for the random background.")
    args = p.parse_args(argv)

    try:
        enzymes = resolve_enzymes(args.enzymes)
    except EnzymeSpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    try:
        ref = pysam.FastaFile(args.reference)
    except (OSError, ValueError) as exc:
        print(f"error: cannot open {args.reference}: {exc}\n"
              f"An .fai index is required (samtools faidx ref.fa).",
              file=sys.stderr)
        return 2
    rng = random.Random(args.seed)
    lengths = dict(zip(ref.references, ref.lengths))

    hits = {e.name: 0 for e in enzymes}
    bg = {e.name: 0 for e in enzymes}
    n = 0
    for chrom, pos in junction_boundaries(args.bam, args.mapq, args.min_jump):
        if pos is None or chrom not in lengths:
            continue
        n += 1
        # background: a random position on the same chromosome, so that base
        # composition and the chromosome mix match the junctions exactly
        rpos = rng.randrange(1000, max(1001, lengths[chrom] - 1000))
        for e in enzymes:
            if motif_at(ref, chrom, pos, e, args.tol):
                hits[e.name] += 1
            if motif_at(ref, chrom, rpos, e, args.tol):
                bg[e.name] += 1
        if args.max_junctions and n >= args.max_junctions:
            break

    if not n:
        print("No ligation junction found. Is the BAM name-sorted, aligned, "
              "and does it carry the MI tag?", file=sys.stderr)
        return 1

    print(f"{n:,} junction boundaries "
          f"(MAPQ >= {args.mapq}, jump > {args.min_jump:,} bp)\n")
    print(f"  {'enzyme':<10} {'site':<8} {'at junctions':>13} "
          f"{'at random':>11} {'enrichment':>11}   verdict")
    for e in enzymes:
        h, b = hits[e.name] / n, max(bg[e.name], 1) / n
        ratio = h / b
        verdict = ("yes, main enzyme" if ratio >= 5
                   else "yes, secondary" if ratio >= 2
                   else "unclear" if ratio >= 1.3
                   else "no")
        print(f"  {e.name:<10} {e.site:<8} {h * 100:>12.1f}% "
              f"{bg[e.name] / n * 100:>10.1f}% {ratio:>10.1f}x   {verdict}")
    print("""
  Read the enrichment column, not the absolute percentages. An enrichment near
  1 means the site is no more common at junctions than anywhere else, so the
  enzyme did not create them.

  The absolute percentages understate the truth, because aligners soft-clip a
  few bases at monomer ends and the boundary drifts off the true cut. Raising
  --tol recovers them (DpnII goes from ~40% to ~79% at --tol 10) while lowering
  every enrichment, since the random background grows too. The ranking between
  enzymes is stable; the thresholds behind the verdict column are a convenience,
  so check the numbers when one lands near a boundary.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
