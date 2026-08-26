# wf-pore-c_AQB — multi-enzyme digestion for Pore-C / CiFi

`pore-c-py digest` accepts **one** restriction enzyme. Many 3C, Hi-C and
Pore-C protocols use **two or more** — DpnII + NlaIII, or the Arima kit's
DpnII + HinfI. Digesting *in silico* with only one of them silently loses
contacts.

This package provides `pore-c-aqb digest`, a drop-in replacement that takes a
list of enzymes and cuts at the union of their sites.

```bash
pore-c-aqb digest DpnII        reads.bam --output monomers.bam   # as before
pore-c-aqb digest DpnII,NlaIII reads.bam --output monomers.bam   # new
```

---

## Why it matters

A concatemer is digested **once**, by every enzyme present in the tube at the
same time. Its fragments are therefore delimited by the union of all
recognition sites.

If the workflow is told about only one enzyme, a monomer that internally spans
a junction made by the *other* enzyme is not split. It contains two distinct
genomic loci, aligns to only one of them, and the remainder is soft-clipped and
discarded. **The contact is lost, silently, with no warning anywhere.**

On the dataset that motivated this tool, junction analysis showed 78–86% of
ligation junctions at `GATC` (DpnII) and 21–26% at `CATG` (NlaIII), while the
pipeline had been run with `--cutter DpnII` alone. Roughly a fifth of the
junctions were invisible to the digest.

## What it reports

The per-enzyme cut report answers a question that is otherwise very hard to
ask: **did this enzyme actually cut?**

```
Digested 4,227 concatemers into 160,688 monomers (156,461 cut points).
  DpnII        GATC             61,229 cuts ( 39.1%)
  NlaIII       CATG             95,424 cuts ( 61.0%)
  (shared)                         192 positions cut by more than one enzyme, counted once
```

An enzyme that failed in the reaction shows up immediately as `0 cuts`, instead
of as a diffuse loss of contacts three analysis steps later. `--stats out.tsv`
writes the same numbers as a table.

---

## Install

```bash
pip install git+https://github.com/AQB/wf-pore-c_AQB.git
```

Requires Python ≥ 3.8, `pysam` and `biopython`. For development:

```bash
git clone https://github.com/AQB/wf-pore-c_AQB.git
cd wf-pore-c_AQB
pip install -e ".[dev]"
pytest
```

## Usage

```
pore-c-aqb digest ENZYME [INPUT] [options]

  ENZYME            One or more enzyme names. Separate with commas:
                    'DpnII', 'DpnII,NlaIII', 'DpnII+NlaIII' all work.
                    Names are Biopython's (Bio.Restriction).
  INPUT             Unaligned BAM of concatemers, '-' for stdin.

  --output PATH     Output BAM ('-' for stdout).
  --header PATH     BAM whose header to copy. Required when reading stdin.
  --stats PATH      Write the per-enzyme cut report as TSV.
  --dry-run         Resolve and print the enzymes, then exit.
  --remove_tags ... Extra SAM tags to strip.
  --max_reads N     Stop after N concatemers.
  --threads N       Threads for BAM compression.
  --debug/--quiet   Logging verbosity.
```

Every option of `pore-c-py digest` is accepted with the same meaning, so an
existing command line keeps working.

## Wiring it into wf-pore-c

See [`docs/INTEGRATION.md`](docs/INTEGRATION.md) and
[`patches/`](patches/). In short, the workflow's digest call becomes:

```diff
- pore-c-py digest "$cutter" ...
+ pore-c-aqb digest "$cutter" ...
```

and `--cutter` may then be given a comma-separated list.

---

## Correctness

Two claims are tested, not asserted.

**1. With one enzyme, output is identical to `pore-c-py`.**

`tests/test_equivalence.py` runs the published `pore-c-py` 2.0.6 inside its own
Docker image and diffs the resulting BAMs record by record — names, sequences,
qualities, `MI`, `Xc`, `MM` and `ML` tags.

Validated on **4,227 real PacBio 3C reads**: 65,456 monomers on both sides,
**zero differing records**, including the 51,332 monomers carrying methylation
tags. The MM/ML recomputation is vendored verbatim from upstream
(`src/pore_c_aqb/_vendored.py`) precisely so that this holds.

**2. With several enzymes, the digest is a strict refinement.**

Every boundary of the single-enzyme digest survives, monomers tile each read
exactly, and the reads reconstruct byte-for-byte from their monomers. Checked
on real data (0 divergences over 4,227 reads) and as a property test over
random sequences.

### A reproducibility bug found along the way

`Bio.Restriction.search()` changed its behaviour at sequence boundaries:

```python
NlaIII.search(Seq("CATGCCTTCTGTGCGAGCCC"))
# biopython 1.82 -> [5]     (the version inside wf-pore-c)
# biopython 1.88 -> []
```

Biopython 1.88 drops sites whose second-strand cut lands exactly on a sequence
edge. About 0.6% of enzyme/sequence pairs are affected — always at position 0
or at the very end, always producing a few-base monomer that never aligns.

The scientific impact is nil, and on a 4,227-read real library the two versions
happened to agree exactly. But the divergence is real and reproducible on
constructed input, which means **the same BAM can give different monomer counts
on two machines** — not acceptable in a tool meant to be cited. This package is
verified identical across both versions on the same real library
(160,688 monomers either way).

`src/pore_c_aqb/sites.py` therefore locates sites itself for palindromic
enzymes — every enzyme used in 3C/Hi-C — under one documented rule, so the
result no longer depends on the installed Biopython. Verified over **2,400
comparisons across 20 palindromic enzymes: zero disagreements outside the
documented boundary case**.

Non-palindromic Type IIS enzymes (BsrI, AcuI, …) cut at a *distance* from their
site, the reverse-strand geometry is fiddly, and none is used to build a Hi-C
library. Those are handed back to Biopython rather than reimplemented on
under-tested arithmetic. The version-dependent edge case can still occur for
them; that is documented rather than hidden.

## Performance

Measured on 4,227 real 3C reads (mean 4.4 kb), single core:

| | time | throughput |
|---|---:|---:|
| `pore-c-py digest DpnII` (upstream) | ~10.8 s | 390 reads/s |
| `pore-c-aqb digest DpnII` | ~10.0 s | 430 reads/s |
| `pore-c-aqb digest DpnII,NlaIII` | ~15.0 s | 282 reads/s |

No regression against upstream. A second enzyme costs about 50% more, split
between the extra site search and the 2.5× larger output.

A single-pass regular expression over all motifs was prototyped and
**rejected**: it gave 243 wrong answers out of 1,600 because an alternation
scan skips overlapping sites (in `CATGATC` the `CATG` match consumes the `GATC`
that starts at offset 3), for only 1.2–1.4× on the realistic multi-enzyme case.
Correctness won. The digest is not the pipeline bottleneck anyway — alignment
dominates by an order of magnitude.

## Design notes

Longer discussion of the choices, with the measurements behind them, in
[`docs/DESIGN.md`](docs/DESIGN.md).

---

## Licence and attribution

Derived from [`pore-c-py`](https://github.com/epi2me-labs/pore-c-py) 2.0.6 by
Oxford Nanopore Technologies PLC, and distributed under the same **Oxford
Nanopore Technologies PLC Public License Version 1.0** (see `LICENSE`).

That licence restricts use to **Research Purposes**. It is not an OSI-approved
open-source licence; read it before depending on this in a commercial setting.

The list of files taken from upstream, and every change made to them, is in
[`NOTICE`](NOTICE).
