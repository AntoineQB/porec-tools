# Changelog

## 1.0.0

First release. Derived from pore-c-py 2.0.6.

### Added
- `pore-c-aqb digest` accepts several restriction enzymes, comma-separated.
  Cut points are the union of all their sites.
- Per-enzyme site report on stderr and via `--stats FILE`: how many recognition
  sites each enzyme has in the reads, its share of the cut points, and the
  observed site spacing next to the spacing expected by chance.
- Enzyme table printed before the digest, and on its own with `--dry-run`:
  site, cut position drawn on both strands (`N^GATC_N`), sticky end, expected
  site spacing. Enzymes leaving incompatible ends are flagged as such.
- `pore-c-aqb enzymes`: lists the enzymes you can digest with, searchable by
  name or by recognition site, with the cut drawn on both strands and the
  sticky end each one leaves. A curated 3C/Hi-C shortlist by default, `--all`
  for the 729 usable ones.
- `pore-c-aqb merge` (also `pore-c-aqb-merge`): glues back the fragments the in
  silico digest split but the enzyme never cut. Outputs merged fragments as
  TSV, contacts as 4DN `.pairs`, and a before/after breakdown of cis distances
  so the artefact is visible. On the reference library, 2,278,102 aligned
  monomers collapse to 450,213 fragments and cis contacts under 10 kb fall from
  79.8% to 26.4% of all pairs.
- `pore-c-aqb junctions` (also `pore-c-aqb-junctions`): the analysis that
  actually determines which enzymes cut a library, from aligned monomers.
- Upstream options the workflow relies on: `--max_monomers`, `--excluded_list`,
  `--excluded_bam`, `--recursive`, `--glob`, and both positional orders for
  the input and the enzyme (wf-pore-c uses one order in each of its two digest
  branches).
- `@PG` provenance line in the output BAM header.
- Enzyme names validated up front, with suggestions on typos. Enzymes with an
  undefined cut position (334 of them) or that cut twice (25) are rejected with
  a clear message rather than mis-digested.

### Fixed
- Digestion no longer depends on the installed Biopython version. `search()`
  changed its behaviour at sequence boundaries between 1.82 and 1.88, which
  made monomer counts differ between machines for the same input.

### Notes on the report
The per-enzyme count is labelled **sites found**, not "cuts". A pre-release
version called it "cuts" and reported, on a real library, `HindIII 1,497 cuts
(9.2%)` — for an enzyme that demonstrably never cut it. The count is of motif
occurrences in the reads, and every motif occurs in genomic DNA by chance; a
failed enzyme still scores highly. The report now prints the chance rate beside
the observed one and states plainly that it is not a test of cutting, pointing
at `pore-c-aqb-junctions`, which is.

### Guarantees
- With a single enzyme, output is identical to `pore-c-py` 2.0.6 — the version
  inside the wf-pore-c container. Verified against that image on 16,758 real
  PacBio 3C reads: 259,214 monomers, 203,342 of them carrying base-modification
  tags, with **zero differences on every field and every tag, value types
  included**.
- Both of wf-pore-c's digest invocations replay end to end through
  `samtools fastq` and `minimap2`. Adding NlaIII to DpnII cuts soft-clipped
  bases from 218,363 to 8,085 on the same reads.
