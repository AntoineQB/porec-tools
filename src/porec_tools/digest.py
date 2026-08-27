# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Multi-enzyme digestion of unaligned concatemers.

Derived from ``pore_c_py.digest`` (Oxford Nanopore Technologies PLC).

The only behavioural change is that :func:`find_cut_points` takes a *list* of
enzymes and returns the sorted union of their cut positions. Everything
downstream (interval splitting, sequence and quality trimming, MM/ML tag
recomputation, monomer naming and tagging) is unchanged, which is what lets
:mod:`tests.test_equivalence` assert byte-identical output against the original
in the single-enzyme case.

Why the union, and not one pass per enzyme
------------------------------------------
A concatemer is digested once, by all the enzymes present in the reaction at
the same time. The resulting fragments are therefore delimited by the union of
all recognition sites. Digesting separately and merging afterwards would be
wrong: it would produce overlapping monomer sets rather than one partition.
"""
from __future__ import annotations

import copy
import itertools
import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterator, Sequence

from Bio.Seq import Seq

from porec_tools import _vendored
from porec_tools.enzymes import ResolvedEnzyme
from porec_tools.sites import find_cuts_for_enzyme

__all__ = ["DigestStats", "find_cut_points", "digest_sequence",
           "digest_sequence_tagged", "get_concatemer_seqs"]

logger = logging.getLogger("porec")


@dataclass
class DigestStats:
    """Counters gathered while digesting.

    ``sites_per_enzyme`` counts recognition sites found in the reads. It
    confirms each enzyme was applied and shows its share of the fragmentation.
    It is deliberately not called "cuts": a site occurring in the sequence is
    not evidence the enzyme cut it, since every motif occurs in genomic DNA by
    chance. See :mod:`porec_tools.report`.
    """

    n_concatemers: int = 0
    n_monomers: int = 0
    #: concatemers dropped for having more than --max_monomers monomers
    n_excluded: int = 0
    n_cut_points: int = 0
    #: total bases read, used to express site density as "1 per N bp"
    n_bases: int = 0
    sites_per_enzyme: Counter = field(default_factory=Counter)
    #: positions carrying a site for more than one enzyme, cut once
    n_shared_cut_points: int = 0

    def summary_lines(self, enzymes: Sequence[ResolvedEnzyme]) -> list[str]:
        """Human-readable report; see :mod:`porec_tools.report`."""
        from porec_tools.report import format_report
        return format_report(self, enzymes)


def find_cut_points(
    sequence: str | Seq,
    enzymes: Sequence[ResolvedEnzyme],
    stats: DigestStats | None = None,
) -> list[int]:
    """Sorted, de-duplicated 0-based cut positions for all enzymes.

    Positions follow ``pore-c-py``'s convention: 0-based, the cut falling
    immediately before the returned index. Site location is delegated to
    :func:`porec_tools.sites.find_cuts_for_enzyme`, which is reproducible
    across Biopython versions for palindromic enzymes.

    Two enzymes can cut at the same position. That is one cut in the tube, so
    it must be one cut here, hence the set.
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence)

    cuts: set[int] = set()
    per_enzyme: dict[str, set[int]] = {}
    for enz in enzymes:
        found = find_cuts_for_enzyme(sequence, enz)
        per_enzyme[enz.name] = found
        cuts |= found

    if stats is not None:
        stats.n_bases += len(sequence)
        for name, found in per_enzyme.items():
            stats.sites_per_enzyme[name] += len(found)
        # a position claimed by k enzymes is counted k times above but once
        # in the digest; record the excess so the report stays honest
        total_claims = sum(len(v) for v in per_enzyme.values())
        stats.n_shared_cut_points += total_claims - len(cuts)
        stats.n_cut_points += len(cuts)

    return sorted(cuts)


def digest_sequence(align, enzymes: Sequence[ResolvedEnzyme],
                    tags_remove=None, stats: DigestStats | None = None,
                    max_monomers=None, legacy_mod_tags: bool = False):
    """Monomers of one concatemer.

    A concatemer excluded by ``max_monomers`` yields nothing here; use
    :func:`digest_sequence_tagged` if you need the excluded read itself.
    """
    for is_monomer, read in digest_sequence_tagged(
            align, enzymes, tags_remove, stats, max_monomers, legacy_mod_tags):
        if is_monomer:
            yield read


def digest_sequence_tagged(align, enzymes: Sequence[ResolvedEnzyme],
                           tags_remove=None, stats: DigestStats | None = None,
                           max_monomers=None, legacy_mod_tags: bool = False):
    """Split one concatemer, yielding ``(is_monomer, read)``.

    Mirrors pore-c-py exactly, with two additions: several enzymes, and the
    ``max_monomers`` exclusion that upstream added after 2.0.6. A read with
    more than ``max_monomers`` pieces is almost always a chimera or an
    over-digested artefact, so upstream drops it whole rather than emitting the
    pieces; ``(False, original_read)`` is yielded so the caller can record it.
    """
    # the move tag massively bloats files, and we don't care for
    # it or handle it in trimming, so force its removal by default.
    if tags_remove is None:
        tags_remove = {'mv'}

    concatemer_id = align.query_name
    cut_points = find_cut_points(align.query_sequence, enzymes, stats)
    read_length = len(align.query_sequence)
    num_digits = len(str(read_length))
    intervals = _vendored.splits_to_intervals(cut_points, read_length)
    num_intervals = len(intervals)

    if max_monomers is not None and num_intervals > max_monomers:
        # bail out before building monomers that would only be discarded
        logger.warning(
            "Dropping read %s, has %d monomers.", concatemer_id, num_intervals)
        yield False, copy.copy(align)
        return

    for idx, (start, end) in enumerate(intervals):
        read = copy.copy(align)
        # trim the sequence and quality
        seq = align.query_sequence[start:end]
        qual = None
        if align.query_qualities:
            qual = align.query_qualities[start:end]
        read.query_sequence = seq
        read.query_qualities = qual
        # deal with mods, upgrading tag from interim to approved spec
        if ('Mm' in tags_remove) or ('MM' in tags_remove):
            for tag in ('Mm', 'Ml', 'MM', 'ML', 'MN'):
                read.set_tag(tag, None)
        elif legacy_mod_tags:
            # byte-identical to pore-c-py 2.0.6: ML as a comma separated string
            mm, ml = _vendored.get_subread_modified_bases(align, start, end)
            for tag in ('Mm', 'Ml'):
                read.set_tag(tag, None)
            read.set_tag("MM", mm)
            read.set_tag("ML", ml)
        else:
            mm, ml = _vendored.get_subread_modified_bases_spec(align, start, end)
            for tag in ('Mm', 'Ml', 'MM', 'ML', 'MN'):
                read.set_tag(tag, None)
            if len(ml):
                read.set_tag("MM", mm)
                read.set_tag("ML", ml)
                read.set_tag("MN", len(seq))
        # lexographically sortable monomer ID
        read.query_name = \
            f"{concatemer_id}:{start:0{num_digits}d}:{end:0{num_digits}d}"
        read.set_tag(
            _vendored.MONOMER_DATA_TAG,
            [start, end, read_length, idx, num_intervals])
        _vendored.set_monomer_data(
            read, start, end, read_length, idx, num_intervals)
        read.set_tag(_vendored.CONCATEMER_ID_TAG, concatemer_id, "Z")
        yield True, read


def get_concatemer_seqs(
    input_file,
    enzymes: Sequence[ResolvedEnzyme],
    remove_tags=None,
    stats: DigestStats | None = None,
    max_monomers=None,
    on_excluded=None,
    max_reads=None,
    legacy_mod_tags: bool = False,
) -> Iterator:
    """Digest concatemers into unaligned monomers.

    :param input_file: pysam.AlignmentFile input
    :param enzymes: enzymes resolved by :func:`porec_tools.enzymes.resolve_enzymes`
    :param remove_tags: additional SAM tags to strip
    :param stats: optional :class:`DigestStats` to accumulate counters into
    :param max_monomers: drop concatemers cut into more pieces than this
    :param on_excluded: called with each dropped read, for --excluded_list/bam
    :param max_reads: read only the first N concatemers
    :param legacy_mod_tags: emit ML as a string, as pore-c-py 2.0.6 did
    """
    if stats is None:
        stats = DigestStats()
    tags_remove = {"mv"}
    if remove_tags:
        tags_remove.update(set(remove_tags))
    stream = input_file.fetch(until_eof=True)
    if max_reads is not None:
        # islice so the N+1st concatemer is never consumed or counted, which
        # is what upstream does and what makes --max_reads N mean exactly N
        stream = itertools.islice(stream, max_reads)
    for align in stream:
        stats.n_concatemers += 1
        for is_monomer, read in digest_sequence_tagged(
                align, enzymes, tags_remove, stats, max_monomers,
                legacy_mod_tags):
            if is_monomer:
                stats.n_monomers += 1
                yield read
            else:
                stats.n_excluded += 1
                if on_excluded is not None:
                    on_excluded(read)
