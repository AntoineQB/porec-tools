# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Walking an aligned monomer BAM, concatemer by concatemer.

The digest names each monomer ``<read>:<start>:<end>`` and tags it with ``MI``
(the concatemer it came from) and ``Xc`` (where it sat along that read). Both
:mod:`pore_c_aqb.merge` and :mod:`pore_c_aqb.junctions` need the same thing:
the monomers of one concatemer, in the order they appeared along the read.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pysam

__all__ = ["read_span", "read_order_key", "concatemer_id", "open_bam",
           "iter_concatemers"]

logger = logging.getLogger("pore-c-aqb")

#: what read_span returns when a monomer's place along the read is unknown
UNKNOWN_SPAN = (-1, -1)


def read_span(aln) -> tuple[int, int]:
    """``(start, end)`` of this monomer along its concatemer.

    ``Xc`` is the ``(start, end, read_length, index, count)`` tag written by the
    digest. Falling back to the name suffix keeps things working on BAMs whose
    tags were stripped, since monomer names end in ``:<start>:<end>``.

    Returns :data:`UNKNOWN_SPAN` when neither is available - a BAM from some
    other source. Callers keep the order the records came in, which is the best
    available guess for a name-grouped file.
    """
    try:
        tag = aln.get_tag("Xc")
        return int(tag[0]), int(tag[1])
    except (KeyError, TypeError, IndexError, ValueError):
        pass
    parts = aln.query_name.rsplit(":", 2)
    if len(parts) == 3:
        try:
            return int(parts[1]), int(parts[2])
        except ValueError:
            pass
    return UNKNOWN_SPAN


def read_order_key(aln) -> int:
    """Where this monomer starts along its concatemer, or -1 if unknown."""
    return read_span(aln)[0]


def concatemer_id(aln) -> str:
    """The read a monomer came from."""
    if aln.has_tag("MI"):
        return aln.get_tag("MI")
    return aln.query_name.rsplit(":", 2)[0]


def open_bam(path: str):
    """Open a BAM, failing with a message a user can act on."""
    if not Path(path).exists():
        raise SystemExit(f"Input file not found: {path}")
    try:
        return pysam.AlignmentFile(str(path), "rb", check_sq=False)
    except (OSError, ValueError) as exc:
        raise SystemExit(
            f"Cannot read {path} as a BAM: {exc}\n"
            f"Is it really a BAM? 'samtools quickcheck' will tell you.")


def iter_concatemers(bam_path: str, require_sorted: bool = True,
                     tell: dict | None = None):
    """Yield ``(concatemer_id, [alignments])``, monomers in read order.

    The BAM must be grouped by read name, which is what the workflow produces
    (``*.ns.bam``). Unmapped, secondary and supplementary records are dropped.

    Pass a dict as ``tell`` to receive the file handle's ``tell`` method under
    that key, so a progress bar can report real progress through the file.

    No MAPQ filter happens here, on purpose. Merging contiguous fragments
    has to run on the complete chain: removing a middle monomer first would
    break the contiguity, and two genomically adjacent blocks would then look
    like two distinct loci. Filter after merging, on the merged blocks.
    """
    bam = open_bam(bam_path)
    warned_unknown = False
    if tell is not None:
        # let a caller watch how far into the file we are, for a progress bar
        tell["tell"] = bam.tell
    try:
        if require_sorted:
            hd = bam.header.to_dict().get("HD", {})
            if hd.get("SO") == "coordinate" and hd.get("GO") != "query":
                raise SystemExit(
                    f"{bam_path} is coordinate-sorted. Monomers of one "
                    f"concatemer must be adjacent; use the name-sorted BAM "
                    f"(the workflow's *.ns.bam) or run "
                    f"'samtools sort -n' first.")
        current = None
        buffer: list = []
        seen: set = set()

        def ordered(items):
            """Monomers in read order; file order where that is unknown."""
            return [a for _, _, a in sorted(
                items, key=lambda t: (t[0] if t[0] >= 0 else t[1], t[1]))]

        for index, aln in enumerate(bam):
            mi = concatemer_id(aln)
            if mi != current:
                if current is not None:
                    yield current, ordered(buffer)
                if mi in seen:
                    raise SystemExit(
                        f"Concatemer {mi} appears in two separate blocks: the "
                        f"BAM is not grouped by read name. Run "
                        f"'samtools sort -n' first.")
                seen.add(mi)
                current, buffer = mi, []
            if (aln.is_unmapped or aln.is_secondary
                    or aln.is_supplementary):
                continue
            start = read_order_key(aln)
            if start < 0 and not warned_unknown:
                warned_unknown = True
                logger.warning(
                    "%s: no Xc tag and no ':start:end' in the read names, so "
                    "the position of each monomer along its read is unknown. "
                    "Falling back to the order records appear in the file. "
                    "Is this BAM really from 'pore-c-aqb digest'?", bam_path)
            buffer.append((start, index, aln))
        if current is not None:
            yield current, ordered(buffer)
    finally:
        bam.close()
