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

On the dataset that motivated this tool, `GATC` (DpnII) is 30× enriched at
ligation junctions and `CATG` (NlaIII) 3–4×, while the pipeline had been run
with `--cutter DpnII` alone. A fifth of the junctions were invisible to the
digest. `pore-c-aqb-junctions` is the analysis that establishes this from
aligned data.

## What it reports

**Before reading any data**, `--dry-run` shows what each enzyme will do — both
strands, so the sticky end is visible:

```
$ pore-c-aqb digest DpnII,NlaIII --dry-run
Enzymes resolved from 'DpnII,NlaIII':
  enzyme  site  cut (both strands)  sticky end  1 site per  recognition
  DpnII   GATC  N^GATC_N            5' GATC     256 bp      palindromic
  NlaIII  CATG  _CATG^              3' CATG     256 bp      palindromic
  These enzymes leave different ends, so they do not ligate to each other:
  each junction should carry a single enzyme's motif (unless the protocol
  fills in ends before ligation).
```

**After the digest**, the per-enzyme table shows how much each enzyme
contributed to the fragmentation:

```
Digested 1,036 concatemers (5,250,464 bases) into 40,307 monomers,
cutting at 39,271 distinct positions.

Recognition sites found in the reads:
  enzyme   site    sites found  % of cuts  1 per (observed)  1 per (chance)
  DpnII    GATC         14,721      37.5%            357 bp          256 bp
  NlaIII   CATG         23,104      58.8%            227 bp          256 bp
  HindIII  AAGCTT        1,497       3.8%          3,507 bp        4,096 bp
```

`--stats out.tsv` writes the same numbers as a table.

### What this table does *not* say

It counts **sites present in the reads**, not cuts. Every motif occurs in
genomic DNA by chance — at roughly the rate in the last column — whether or not
the enzyme was ever in the tube. In the run above, HindIII shows 1,497 sites at
1 per 3,507 bp against a chance rate of 1 per 4,096 bp: exactly what you get
from a genome that never saw the enzyme. HindIII did **not** cut this library.

So the count is useful for deciding whether a second enzyme is worth adding,
and for catching a wrong enzyme name (`0 sites` triggers a warning). It is
**not** a test of whether the enzyme worked. That test needs aligned data:
compare the motifs at ligation junctions against random positions on the same
chromosomes, which is what `pore-c-aqb-junctions` does — see
[`docs/INTEGRATION.md`](docs/INTEGRATION.md). A pre-release version of this
tool labelled the column "cuts", which invited precisely the wrong conclusion.

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
pore-c-aqb digest ENZYME [INPUT ...] [options]
pore-c-aqb digest [INPUT ...] ENZYME [options]     # upstream's order

  ENZYME            One or more enzyme names. Separate with commas:
                    'DpnII', 'DpnII,NlaIII', 'DpnII+NlaIII' all work.
                    Names are Biopython's (Bio.Restriction).
  INPUT             Unaligned BAM(s) of concatemers, a directory, or '-'
                    for stdin (the default).

  --output PATH        Output BAM ('-' for stdout).
  --header PATH        BAM whose header to copy. Required when reading stdin.
  --stats PATH         Write the per-enzyme site report as TSV.
  --dry-run            Print what each enzyme will do, then exit.
  --remove_tags ...    Extra SAM tags to strip.
  --max_reads N        Take only the first N concatemers.
  --max_monomers N     Drop a concatemer cut into more than N monomers.
  --excluded_list PATH Names of the reads dropped by --max_monomers.
  --excluded_bam PATH  The dropped reads themselves.
  --recursive          Search an input directory recursively.
  --glob PATTERN       Which files to take from a directory (default *.bam).
  --threads N          Threads for BAM compression.
  --debug/--quiet      Logging verbosity.
```

Both positional orders are accepted because wf-pore-c uses **both**: its
chunked branch runs `digest "$cutter"` on stdin, its other branch runs
`digest concatemers.bam "$cutter"`.

Every option of `pore-c-py digest` is accepted with the same meaning, so an
existing command line keeps working.

A second command ships with the package:

```
pore-c-aqb-junctions ALIGNED_BAM REFERENCE --enzymes DpnII,NlaIII[,...]
```

It reports which enzymes actually cut a library, from aligned monomers. This is
the question the digest report cannot answer; see
[`docs/INTEGRATION.md`](docs/INTEGRATION.md).

## Wiring it into wf-pore-c

See [`docs/INTEGRATION.md`](docs/INTEGRATION.md) and
[`patches/`](patches/). In short, the workflow's digest call becomes:

```diff
- pore-c-py digest "$cutter" ...
+ pore-c-aqb digest "$cutter" ...
```

at **both** of the workflow's digest call sites — it invokes the digest once
per chunking branch, with the positionals in a different order each time — and
`--cutter` may then be given a comma-separated list.

---

## Correctness

Three claims are tested, not asserted.

**1. With one enzyme, output is identical to `pore-c-py`.**

`tests/test_equivalence.py` runs the real `pore-c-py` inside the published
**wf-pore-c** image — the one the workflow actually uses, carrying pore-c-py
**2.0.6** and biopython 1.82 — and diffs the resulting BAMs record by record.
The strict comparison takes every field and **every tag, including each tag's
value type**, so a divergence in something nobody thought to list still fails.

Validated on **16,758 real PacBio 3C reads → 259,214 monomers**, of which
203,342 carry base-modification tags: **zero differences, on every field and
every tag**. The MM/ML recomputation is vendored verbatim from upstream
(`src/pore_c_aqb/_vendored.py`) precisely so that this holds.

> Pin the image, not the package name. The standalone `ontresearch/pore-c-py`
> image is a later 2.1.x which changed the mod-base tags — `ML` became a uint8
> array and `MN` was added. Diffing against it reports differences that are
> upstream's own version bump. This tool tracks 2.0.6 because that is what
> wf-pore-c ships.

**2. With several enzymes, the digest is a strict refinement.**

Every boundary of the single-enzyme digest survives, monomers tile each read
exactly, and the reads reconstruct byte-for-byte from their monomers. Checked
on **16,758 real concatemers**: 0 lost boundaries, 0 tiling gaps, 259,214
monomers becoming 634,479 — and as a property test over random sequences.

**3. The patched workflow actually runs.**

Both of wf-pore-c's digest invocations were replayed on real reads, through
`samtools fastq` and `minimap2`, exactly as the workflow chains them. Adding
the second enzyme does what it is supposed to do:

| `--cutter` | monomers | aligned | soft-clipped bases |
|---|---:|---:|---:|
| `DpnII` | 17,076 | 7,390 | 218,363 |
| `DpnII,NlaIII` | 39,910 | 8,543 | **8,085** |

Soft-clipping falls by 96%. Those clipped bases were the second locus of a
monomer that spanned an undigested NlaIII junction — the contacts the
single-enzyme digest was throwing away.

### A reproducibility bug found along the way

`Bio.Restriction.search()` changed its behaviour at sequence boundaries:

```python
NlaIII.search(Seq("CATGCCTTCTGTGCGAGCCC"))
# biopython 1.82 -> [5]     (the version inside wf-pore-c)
# biopython 1.88 -> []
```

Biopython 1.88 drops sites whose second-strand cut lands exactly on a sequence
edge — always at position 0 or at the very end, always producing a few-base
monomer that never aligns. It is rare (2 cases in 350 pairs constructed to hit
boundaries; 0 in 2,400 random sequences), but it is deterministic, not random:
the same input reproduces it every time.

The scientific impact is nil, and on a real library the two versions happen to
agree exactly. But the divergence is real and reproducible on constructed
input, which means **the same BAM can give different monomer counts on two
machines** — not acceptable in a tool meant to be cited. This package is
verified identical across both versions on the same real library: 16,758 reads,
`DpnII,NlaIII`, **634,479 monomers either way, zero differing records**. The
full test suite (170 tests) runs on both versions in CI.

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

Measured on 16,758 real 3C reads (mean 4.4 kb), single core, second run of two:

| | time | throughput |
|---|---:|---:|
| `pore-c-py digest DpnII` (upstream, in its container) | 42.3 s | 395 reads/s |
| `pore-c-aqb digest DpnII` | 39.9 s | 420 reads/s |
| `pore-c-aqb digest DpnII,NlaIII` | 60.8 s | 275 reads/s |

No regression against upstream. A second enzyme costs about 50% more, split
between the extra site search and the 2.4× larger output (259,214 → 634,479
monomers). The digest is not the pipeline bottleneck: alignment of the same
reads takes an order of magnitude longer.

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
