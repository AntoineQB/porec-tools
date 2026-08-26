# Changelog

## 1.0.0

First release. Derived from pore-c-py 2.0.6.

### Added
- `pore-c-aqb digest` accepts several restriction enzymes, comma-separated.
  Cut points are the union of all their sites.
- Per-enzyme cut report on stderr and via `--stats FILE`, so an enzyme that did
  not cut is visible immediately.
- `--dry-run` to check an enzyme specification without reading any data.
- `@PG` provenance line in the output BAM header.
- Enzyme names validated up front, with suggestions on typos.

### Fixed
- Digestion no longer depends on the installed Biopython version. `search()`
  changed its behaviour at sequence boundaries between 1.82 and 1.88, which
  made monomer counts differ between machines for the same input.

### Guarantees
- With a single enzyme, output is byte-identical to `pore-c-py` 2.0.6.
  Verified against the upstream Docker image on 4,227 real PacBio 3C reads:
  65,456 monomers, zero differing records, methylation tags included.
