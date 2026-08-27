# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Resolution and validation of one or more restriction enzymes.

This module is the heart of the multi-enzyme extension. ``pore-c-py`` accepts a
single enzyme name; many 3C/Hi-C protocols use two or more (for instance
DpnII + NlaIII, or the Arima kit's DpnII + HinfI). Digesting with only one of
them leaves the other's ligation junctions undetected: the affected monomer
spans two distinct loci, aligns to only one of them, and the remainder is
soft-clipped away. The contact is silently lost.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from Bio import Restriction
from Bio.Seq import Seq

__all__ = [
    "COMMON_3C_ENZYMES",
    "iter_usable",
    "EnzymeSpecError",
    "ResolvedEnzyme",
    "parse_enzyme_spec",
    "resolve_enzymes",
    "describe_enzymes",
]


class EnzymeSpecError(ValueError):
    """Raised when an enzyme specification cannot be resolved."""


@dataclass(frozen=True)
class ResolvedEnzyme:
    """A restriction enzyme validated for use in a digest."""

    name: str
    site: str
    fst5: int
    is_palindromic: bool
    enzyme: object = field(repr=False)

    def search(self, seq: Seq) -> list[int]:
        """Cut positions, 1-based, as returned by Biopython."""
        return self.enzyme.search(seq)


def parse_enzyme_spec(spec: str | Sequence[str]) -> list[str]:
    """Split a user-supplied enzyme specification into names.

    Accepts ``"DpnII"``, ``"DpnII,NlaIII"``, ``"DpnII+NlaIII"``,
    ``"DpnII NlaIII"`` or an already-split sequence. Order is preserved and
    duplicates are dropped, because a site cut twice is still one cut.
    """
    if not isinstance(spec, str):
        raw: Iterable[str] = spec
    else:
        normalised = spec.replace("+", ",").replace(";", ",").replace(" ", ",")
        raw = normalised.split(",")

    names: list[str] = []
    for item in raw:
        name = item.strip()
        if not name:
            continue
        if name not in names:
            names.append(name)
    if not names:
        raise EnzymeSpecError(
            "No enzyme given. Provide one or more names, e.g. 'DpnII' or "
            "'DpnII,NlaIII'."
        )
    return names


def _suggest(name: str, limit: int = 4) -> list[str]:
    """Enzyme names close to a mistyped one (case-insensitive first)."""
    known = [n for n in dir(Restriction) if not n.startswith("_")]
    exact_ci = [n for n in known if n.lower() == name.lower()]
    if exact_ci:
        return exact_ci
    return difflib.get_close_matches(name, known, n=limit, cutoff=0.6)


def _resolve_one(name: str) -> ResolvedEnzyme:
    enzyme = getattr(Restriction, name, None)
    if enzyme is None or not hasattr(enzyme, "search"):
        hint = ""
        suggestions = _suggest(name)
        if suggestions:
            hint = f" Did you mean: {', '.join(suggestions)}?"
        raise EnzymeSpecError(f"Unknown restriction enzyme: {name!r}.{hint}")

    # Enzymes that cut on both sides of their recognition site would produce
    # two cut points per site; pore-c-py does not model this and neither do we.
    if enzyme.cut_twice():
        raise EnzymeSpecError(
            f"Enzyme {name!r} cuts twice, which is not supported."
        )
    # A blunt/undefined cutting position cannot be turned into a split point.
    if enzyme.fst5 is None:
        raise EnzymeSpecError(
            f"Enzyme {name!r} has no defined 5' cutting position."
        )

    return ResolvedEnzyme(
        name=name,
        site=str(enzyme.site),
        fst5=int(enzyme.fst5),
        is_palindromic=bool(enzyme.is_palindromic()),
        enzyme=enzyme,
    )


def resolve_enzymes(spec: str | Sequence[str]) -> list[ResolvedEnzyme]:
    """Resolve a specification into validated enzymes.

    Every name is validated before any is used, so a typo in the second enzyme
    fails immediately rather than after hours of digestion.
    """
    names = parse_enzyme_spec(spec)
    resolved = [_resolve_one(n) for n in names]

    # Two names for the same recognition sequence and cut position (DpnII and
    # MboI, say) would do identical work twice.
    seen: dict[tuple[str, int], str] = {}
    unique: list[ResolvedEnzyme] = []
    for enz in resolved:
        key = (enz.site, enz.fst5)
        if key in seen:
            continue
        seen[key] = enz.name
        unique.append(enz)
    return unique


def describe_enzymes(enzymes: Sequence[ResolvedEnzyme]) -> str:
    """One-line human-readable summary, for logs and BAM headers."""
    return ", ".join(f"{e.name}({e.site})" for e in enzymes)


#: Enzymes that actually turn up in 3C / Hi-C / Pore-C protocols. Not a
#: restriction on what you may use - any name below is accepted - just the
#: short list worth showing first, since Biopython knows 729 usable enzymes and
#: a wall of them helps nobody.
COMMON_3C_ENZYMES = [
    "DpnII", "MboI", "Sau3AI",      # GATC, the classic 4-cutter
    "NlaIII",                       # CATG, the usual second enzyme
    "HinfI",                        # GANTC, the Arima kit's second enzyme
    "MseI", "MluCI", "CviQI", "BfaI", "AluI", "MspI", "TaqI", "DdeI",
    "HindIII", "BglII", "EcoRI", "BamHI", "NcoI", "KpnI", "XhoI",
    "NotI",                         # 8-cutter, for very coarse maps
]


def iter_usable(query: str | None = None, common_only: bool = False):
    """Enzymes this tool can digest with, optionally filtered.

    ``query`` matches a name (case-insensitively, as a substring) or a
    recognition site. Enzymes Biopython cannot place a cut for, and those that
    cut twice, are left out: they are rejected by :func:`resolve_enzymes`
    anyway, so listing them would only invite a confusing error.
    """
    names = COMMON_3C_ENZYMES if common_only else sorted(
        Restriction.AllEnzymes.as_string(), key=str.lower)
    needle = query.upper() if query else None
    for name in names:
        try:
            enzyme = _resolve_one(name)
        except EnzymeSpecError:
            continue
        if needle and needle not in name.upper() and needle != enzyme.site:
            continue
        yield enzyme
