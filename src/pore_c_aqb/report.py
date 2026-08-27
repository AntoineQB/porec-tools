# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Formatting of the enzyme table and the digest report.

What the report can and cannot tell you
---------------------------------------
It counts **recognition sites found in the reads**, per enzyme. That is a
factual, useful number: it confirms the enzyme was really applied, catches a
name typo, and shows how much each enzyme contributes to the fragmentation.

It is *not* evidence that the enzyme cut in the tube. Every enzyme's site
occurs in genomic DNA whether or not the enzyme was in the reaction — AAGCTT
turns up roughly every 4 kb by chance alone. An enzyme that failed completely
still shows a large site count here. The report therefore says "sites found",
never "cuts", and prints the chance rate beside the observed one so the reader
can see that the two are expected to match. What the numbers do tell you is how
much fragmentation each enzyme can possibly contribute, which is what you need
when deciding whether adding an enzyme was worth it.

Telling a working enzyme from a dud requires the **aligned** data: take the
junctions between consecutive monomers that jump in the genome, and test which
motif sits at those boundaries against a random-position background. See
docs/INTEGRATION.md.

Displaying the cut
------------------
``fst5=0`` means nothing to a bench biologist. ``N^GATC_N`` does: it shows both
strands, so the sticky end is visible. That string comes from Biopython's
``elucidate()``, which was verified identical between biopython 1.82 and 1.88
across all 1,086 shared enzymes — unlike ``search()``, which is not (see
``sites.py``). It is used for display only; cut positions never come from it.
"""
from __future__ import annotations

from typing import Sequence

from pore_c_aqb.enzymes import ResolvedEnzyme

__all__ = [
    "expected_site_spacing",
    "describe_cut",
    "describe_overhang",
    "enzyme_table",
    "format_report",
    "list_enzymes",
]


def expected_site_spacing(site: str) -> float:
    """Mean distance between occurrences of a site in random DNA, in bases.

    A fully specified base divides the expected spacing by 4; an ambiguity code
    divides it by 4 / (number of bases it allows), so GANTC gives 256 rather
    than 1024 — the N is free.
    """
    from pore_c_aqb.sites import IUPAC_CLASSES

    spacing = 1.0
    for char in site.upper():
        cls = IUPAC_CLASSES.get(char, char)
        n_bases = sum(cls.count(b) for b in "ACGT") if len(cls) > 1 else 1
        spacing *= 4.0 / n_bases
    return spacing


def describe_cut(enzyme: ResolvedEnzyme) -> str:
    """Both strands around the cut, e.g. ``N^GATC_N``.

    ``^`` is the top-strand cut, ``_`` the bottom-strand one. Falls back to a
    single-strand form if Biopython cannot elucidate the enzyme; enzymes it
    cannot handle at all are rejected upstream by ``resolve_enzymes``.
    """
    try:
        text = enzyme.enzyme.elucidate()
        if "^" in text and "sorry" not in text:
            return text
    except Exception:                    # pragma: no cover - defensive
        pass
    fst5 = enzyme.fst5
    if fst5 is not None and 0 <= fst5 <= len(enzyme.site):
        return f"{enzyme.site[:fst5]}^{enzyme.site[fst5:]}"
    return f"{enzyme.site} (cuts at offset {fst5})"


def describe_overhang(enzyme: ResolvedEnzyme) -> str:
    """The single-stranded end left behind, e.g. ``5' GATC`` or ``blunt``.

    Two enzymes leaving different overhang types cannot ligate to each other,
    which is why a junction in such a library carries one enzyme's motif and
    not a hybrid. This holds for direct sticky-end ligation; a protocol that
    fills in the ends before ligation makes everything blunt and compatible.
    """
    ovhg = getattr(enzyme.enzyme, "ovhg", None)
    if ovhg is None:
        return "?"
    if ovhg == 0:
        return "blunt"
    seq = getattr(enzyme.enzyme, "ovhgseq", "") or ""
    end = "5'" if ovhg < 0 else "3'"
    return f"{end} {seq}" if seq else f"{end} {abs(ovhg)} nt"


def enzyme_table(enzymes: Sequence[ResolvedEnzyme],
                 note_compatibility: bool = True) -> list[str]:
    """What each enzyme will do, before any data is read.

    ``note_compatibility`` adds a line when the chosen enzymes leave ends that
    cannot ligate to each other. That is worth saying when the user has picked
    a set to digest with, and meaningless in a catalogue listing, where the
    enzymes shown are alternatives rather than a combination.
    """
    rows = [
        (e.name, e.site, describe_cut(e), describe_overhang(e),
         f"{expected_site_spacing(e.site):,.0f} bp",
         "palindromic" if e.is_palindromic else "non-palindromic")
        for e in enzymes
    ]
    head = ("enzyme", "site", "cut (both strands)", "sticky end",
            "1 site per", "recognition")
    widths = [max(len(h), *(len(r[i]) for r in rows))
              for i, h in enumerate(head)]
    lines = ["  " + "  ".join(h.ljust(w) for h, w in zip(head, widths)).rstrip()]
    for r in rows:
        lines.append(
            "  " + "  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip())

    ends = {describe_overhang(e) for e in enzymes}
    if (note_compatibility and len(enzymes) > 1
            and len(ends) == len(enzymes) and "?" not in ends):
        lines.append(
            "  These enzymes leave different ends, so they do not ligate to "
            "each other: each junction should carry a single enzyme's motif "
            "(unless the protocol fills in ends before ligation).")
    return lines


def format_report(stats, enzymes: Sequence[ResolvedEnzyme]) -> list[str]:
    """The post-digest summary, as lines ready to log."""
    lines = [
        f"Digested {stats.n_concatemers:,} concatemers "
        f"({stats.n_bases:,} bases) into {stats.n_monomers:,} monomers, "
        f"cutting at {stats.n_cut_points:,} distinct positions."
    ]
    if stats.n_excluded:
        lines.append(
            f"Excluded {stats.n_excluded:,} concatemers "
            f"({100.0 * stats.n_excluded / stats.n_concatemers:.1f}%) for "
            f"exceeding --max_monomers.")
    if not stats.n_concatemers or not stats.n_cut_points:
        return lines

    rows = []
    for enz in enzymes:
        n = stats.sites_per_enzyme.get(enz.name, 0)
        rows.append((
            enz.name,
            enz.site,
            f"{n:,}",
            f"{100.0 * n / stats.n_cut_points:.1f}%",
            f"{stats.n_bases / n:,.0f} bp" if n else "never",
            f"{expected_site_spacing(enz.site):,.0f} bp",
        ))
    head = ("enzyme", "site", "sites found", "% of cuts",
            "1 per (observed)", "1 per (chance)")
    widths = [max(len(h), *(len(r[i]) for r in rows))
              for i, h in enumerate(head)]
    align = (str.ljust, str.ljust, str.rjust, str.rjust, str.rjust, str.rjust)

    lines.append("")
    lines.append("Recognition sites found in the reads:")
    lines.append("  " + "  ".join(
        f(h, w) for f, h, w in zip(align, head, widths)).rstrip())
    for r in rows:
        lines.append("  " + "  ".join(
            f(c, w) for f, c, w in zip(align, r, widths)).rstrip())

    if stats.n_shared_cut_points:
        lines.append(
            f"  {stats.n_shared_cut_points:,} positions carry a site for more "
            f"than one enzyme; each is one cut, so the percentages above sum "
            f"to slightly more than 100.")

    silent = [e.name for e in enzymes
              if stats.sites_per_enzyme.get(e.name, 0) == 0]
    if silent:
        lines.append("")
        lines.append(
            f"  WARNING: not a single site for {', '.join(silent)}. The enzyme "
            f"name is probably wrong, or these reads are not from this run.")

    lines.append("")
    lines.append(
        "  These are sites PRESENT in the reads, not proof an enzyme cut. "
        "Every motif occurs in genomic DNA by chance, at about the rate in "
        "the last column, whether or not the enzyme was in the tube - so a "
        "large count here means nothing on its own. To find out which enzyme "
        "really cut, compare motifs at ligation junctions against random "
        "positions in the ALIGNED data; see docs/INTEGRATION.md.")
    return lines


def list_enzymes(query: str | None = None, common_only: bool = True) -> list[str]:
    """The catalogue, as lines ready to print."""
    from pore_c_aqb.enzymes import iter_usable

    found = list(iter_usable(query, common_only=common_only))
    if not found:
        return [f"No usable enzyme matches {query!r}.",
                "Try a name fragment ('Dpn'), a recognition site ('GATC'), "
                "or --all to see every one."]
    lines = enzyme_table(found, note_compatibility=False)
    scope = "commonly used in 3C/Hi-C" if common_only else "usable"
    filtered = f" matching {query!r}" if query else ""
    lines.append("")
    lines.append(f"  {len(found)} enzyme(s) {scope}{filtered}.")
    if common_only and not query:
        lines.append("  This is a short list, not a restriction: any of the "
                     "729 enzymes Biopython can place a cut for is accepted. "
                     "Use --all to see them, or pass a name or site to search.")
    return lines
