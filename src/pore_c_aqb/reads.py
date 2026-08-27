# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Walking an aligned monomer BAM, concatemer by concatemer.

The digest names each monomer ``<read>:<start>:<end>`` and tags it with ``MI``
(the concatemer it came from) and ``Xc`` (where it sat along that read). Both
:mod:`pore_c_aqb.merge` and :mod:`pore_c_aqb.junctions` need the same thing:
the monomers of one concatemer, in the order they appeared along the read.
"""
from __future__ import annotations

import pysam

__all__ = ["read_order_key", "concatemer_id", "iter_concatemers"]


def read_order_key(aln) -> int:
    """Where this monomer starts along its concatemer.

    ``Xc`` is the ``(start, end, read_length, index, count)`` tag written by the
    digest. Falling back to the name suffix keeps things working on BAMs whose
    tags were stripped, since monomer names end in ``:<start>:<end>``.
    """
    try:
        return int(aln.get_tag("Xc")[0])
    except (KeyError, TypeError, IndexError):
        parts = aln.query_name.rsplit(":", 2)
        return int(parts[1]) if len(parts) == 3 else 0


def concatemer_id(aln) -> str:
    """The read a monomer came from."""
    if aln.has_tag("MI"):
        return aln.get_tag("MI")
    return aln.query_name.rsplit(":", 2)[0]


def iter_concatemers(bam_path: str, require_sorted: bool = True):
    """Yield ``(concatemer_id, [alignments])``, monomers in read order.

    The BAM must be grouped by read name, which is what the workflow produces
    (``*.ns.bam``). Unmapped, secondary and supplementary records are dropped.

    **No MAPQ filter happens here, deliberately.** Merging contiguous fragments
    has to run on the complete chain: removing a middle monomer first would
    break the contiguity, and two genomically adjacent blocks would then look
    like two distinct loci. Filter after merging, on the merged blocks.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
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
        for aln in bam:
            mi = concatemer_id(aln)
            if mi != current:
                if current is not None:
                    yield current, sorted(buffer, key=read_order_key)
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
            buffer.append(aln)
        if current is not None:
            yield current, sorted(buffer, key=read_order_key)
    finally:
        bam.close()
