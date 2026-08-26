# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Functions taken unchanged from ``pore-c-py`` 2.0.6.

Reproduced verbatim from ``pore_c_py/digest.py`` and ``pore_c_py/utils.py``
(Oxford Nanopore Technologies PLC, ONT Public License v1.0) so that the output
of this tool is byte-for-byte identical to the original when a single enzyme is
used. Do not "improve" anything here: any divergence breaks that guarantee, and
``tests/test_equivalence.py`` will fail.

Modifications: none. Only the surrounding module layout differs.
"""
from __future__ import annotations

__all__ = [
    "CONCATEMER_ID_TAG",
    "MONOMER_DATA_TAG",
    "get_subread_modified_bases",
    "splits_to_intervals",
]

CONCATEMER_ID_TAG = "MI"
MONOMER_DATA_TAG = "Xc"
WALK_TAG = "Xw"


def get_subread_modified_bases(align, start, end):
    """Get modified bases subread.

    :param align: pysam.AlignedSegment.
    :param start: start coordinate to trim to.
    :param end: exclusive end co-ordinate.
    """
    mm_str, ml_str = "", ""
    base_indices = {}
    seq = align.query_sequence[start:end]
    for mod_key, mod_data in align.modified_bases.items():
        # find the modifications that overlap the subread
        idx = [
            x for x in range(len(mod_data)) if
            start <= mod_data[x][0] < end]
        probs_dic = dict(mod_data)
        if not idx:  # no mod bases (of this type) in the subread
            continue
        try:
            canonical_base, strand, skip_scheme, mod_type = mod_key
        except ValueError:
            canonical_base, strand, mod_type = mod_key
            skip_scheme = ""
        if canonical_base == "N":
            base_indices[canonical_base] = list(range(len(seq)))
        elif canonical_base not in base_indices:
            base_indices[canonical_base] = [
                x for x, b in enumerate(seq) if b.upper() == canonical_base
            ]
        base_offsets, probs = zip(*[mod_data[i] for i in idx])
        strand = "+" if strand == 0 else "-"
        deltas = []
        counter = 0
        probs = []
        for seq_idx in base_indices[canonical_base]:
            orig_idx = seq_idx + start
            if orig_idx in base_offsets:  # is modified
                probs += [probs_dic[orig_idx]]
                deltas.append(str(counter))
                counter = 0
            else:
                counter += 1
            deltas_formatted = ','.join(deltas)
            prob_str = ",".join(map(str, probs))
        mm_str += (
                f"{canonical_base}{strand}{mod_type}{skip_scheme}"
                f",{deltas_formatted};"
            )
        ml_str += f"{prob_str},"
    return mm_str, ml_str


def splits_to_intervals(positions, length):
    """Split to intervals."""
    if len(positions) == 0:
        return [(0, length)]
    prefix, suffix = [], []
    if positions[0] != 0:
        prefix = [0]
    if positions[-1] != length:
        suffix = [length]
    breaks = prefix + positions + suffix
    return [(start, end) for start, end in zip(breaks[:-1], breaks[1:])]


def set_monomer_data(align, start, end, read_length, idx, num_intervals):
    """Set the monomer data on an alignment."""
    align.set_tag(
        MONOMER_DATA_TAG,
        [start, end, read_length, idx, num_intervals])
