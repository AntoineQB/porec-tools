# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Location of restriction cut sites, reproducibly.

Why this module exists
----------------------
``pore-c-py`` calls ``Bio.Restriction.<Enzyme>.search()``. That is correct, but
its behaviour at the *edges* of a sequence changed between Biopython releases::

    NlaIII.search(Seq("CATGCCTTCTGTGCGAGCCC"))
        biopython 1.82  ->  [5]        (the version inside wf-pore-c)
        biopython 1.88  ->  []

Biopython 1.88 drops sites whose second-strand cut would land exactly on a
sequence boundary; 1.82 keeps them. Measured on 350 enzyme/sequence pairs the
two versions disagree on 0.6% of cases, always at position 0 or at the very
end, and the monomers involved are a few bases long and never align. The
scientific impact is nil — but *the same BAM digested on two machines can give
different monomer counts*, which is not acceptable in a tool meant to be
shared and cited.

What we do about it
-------------------
For **palindromic** enzymes — which is every enzyme used in 3C/Hi-C protocols
(DpnII, NlaIII, HinfI, MseI, HindIII, EcoRI, BglII, NcoI, DdeI, Csp6I, ApoI,
BamHI, AluI, HaeIII, TaqI, …) — we locate sites ourselves and apply one
documented rule, so the result no longer depends on the installed Biopython.

For **non-palindromic** enzymes we hand the search back to Biopython. Those are
Type IIS enzymes that cut at a *distance* from their recognition site (AcuI
cuts 16 bases downstream), the reverse-strand geometry is fiddly, and none of
them is used to make a Hi-C library. Shipping under-tested arithmetic there
would be worse than delegating to a battle-tested implementation. The
version-dependent edge case can then still occur for those enzymes; that is
documented rather than hidden.

The rule, for palindromic enzymes
---------------------------------
A cut is emitted at ``site_start + fst5`` for every occurrence of the
recognition site, **including overlapping occurrences**. Cuts outside
``[0, len(seq)]`` are discarded; cuts at 0 or ``len(seq)`` are kept, because
``splits_to_intervals`` already handles them without creating empty monomers.

Overlaps matter: in ``CATGATC`` the NlaIII site starts at 0 and the DpnII site
at 3. A single non-overlapping scan finds one and silently loses the other.
"""
from __future__ import annotations

import re
from functools import lru_cache

from Bio.Seq import Seq

__all__ = ["IUPAC_CLASSES", "site_regex", "find_cuts_for_enzyme"]

#: IUPAC ambiguity codes, as regular-expression character classes.
IUPAC_CLASSES = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def _pattern(site: str) -> str:
    try:
        return "".join(IUPAC_CLASSES[c] for c in site.upper())
    except KeyError as exc:  # pragma: no cover - guarded by enzyme validation
        raise ValueError(
            f"Unsupported character {exc.args[0]!r} in recognition site "
            f"{site!r}"
        ) from None


@lru_cache(maxsize=256)
def site_regex(site: str) -> re.Pattern:
    """Compiled regex matching a recognition site, overlaps included.

    The lookahead makes each match zero-width, so ``finditer`` advances one
    base at a time and reports overlapping occurrences. Compiled once per site
    and reused for every read.
    """
    return re.compile(f"(?=({_pattern(site)}))")


def find_cuts_for_enzyme(sequence: str, enzyme) -> set[int]:
    """0-based cut positions of one enzyme on one sequence.

    :param sequence: the read; characters other than ACGT simply never match.
    :param enzyme: a :class:`pore_c_aqb.enzymes.ResolvedEnzyme`.
    """
    n = len(sequence)

    if not enzyme.is_palindromic:
        # Type IIS and friends: defer to Biopython, see the module docstring.
        return {x - 1 for x in enzyme.search(Seq(sequence))}

    cuts: set[int] = set()
    for match in site_regex(enzyme.site).finditer(sequence.upper()):
        pos = match.start() + enzyme.fst5
        if 0 <= pos <= n:
            cuts.add(pos)
    return cuts
