# wf-pore-c_AQB

**A modified version of Oxford Nanopore's [`pore-c-py`](https://github.com/epi2me-labs/pore-c-py) / [`wf-pore-c`](https://github.com/epi2me-labs/wf-pore-c), for Pore-C and PacBio CiFi libraries.**

This is not a new pipeline. It is the upstream digest, kept byte-for-byte
identical where it was already right, with three things fixed that were
silently costing contacts:

| | Problem | Fix |
|---|---|---|
| 1 | `--cutter` accepts **one** enzyme. Many protocols use two or more, and the second enzyme's junctions are never cut. | `pore-c-aqb digest DpnII,NlaIII` |
| 2 | The *in silico* digest cuts at **every** site, including the many the enzyme never cut. This inflates the Hi-C diagonal with contacts that never happened. | `pore-c-aqb merge` |
| 3 | Nothing tells you which enzymes actually cut your library. | `pore-c-aqb junctions` |

Plus one reproducibility bug in a dependency, described [below](#a-reproducibility-bug-found-along-the-way).

With a single enzyme, the digest output is **identical to upstream** — verified
record by record, tag by tag, against the real tool. See [Correctness](#correctness).

---

## Where it fits

```
   concatemers.bam  (unaligned reads from the sequencer)
          |
          v
   +-----------------------+
   |  pore-c-aqb digest    |   cut each read into monomers
   |  DpnII,NlaIII         |   <- FIX 1: several enzymes
   +-----------------------+
          |  monomers.bam
          v
   +-----------------------+
   |  minimap2             |   each monomer gets a genomic position
   +-----------------------+
          |  aligned.ns.bam
          +------------------------------+
          v                              v
   +-----------------------+   +-----------------------+
   |  pore-c-aqb merge     |   | pore-c-aqb junctions  |
   |  undo the false cuts  |   | which enzymes cut?    |
   |  <- FIX 2             |   | <- FIX 3 (diagnostic) |
   +-----------------------+   +-----------------------+
          |  fragments.tsv.gz + contacts.pairs
          v
      cooler / juicer
```

---

## Install

```bash
pip install git+https://github.com/YOUR-USERNAME/wf-pore-c_AQB.git
```

> **Before the first push**, point the URLs at your own account:
> ```bash
> grep -rl YOUR-USERNAME . --exclude-dir=.git \
>   | xargs sed -i 's/YOUR-USERNAME/your-github-handle/g'
> ```

Needs Python >= 3.8, `pysam` and `biopython` — both pulled in automatically.
Three commands are installed:

```bash
pore-c-aqb enzymes      # which enzymes can I use?
pore-c-aqb digest
pore-c-aqb merge        # also available as pore-c-aqb-merge
pore-c-aqb junctions    # also available as pore-c-aqb-junctions
```

For development:

```bash
git clone https://github.com/YOUR-USERNAME/wf-pore-c_AQB.git
cd wf-pore-c_AQB
pip install -e ".[dev]"
pytest
```

---

## Quick start

```bash
# 0. find your enzymes, then check what the digest will do with them
pore-c-aqb enzymes                    # the 3C/Hi-C shortlist
pore-c-aqb enzymes CATG               # everything cutting CATG
pore-c-aqb digest DpnII,NlaIII --dry-run

# 1. digest with both enzymes
pore-c-aqb digest DpnII,NlaIII reads.bam \
    --output monomers.bam --stats digest_stats.tsv --threads 8

# 2. align (unchanged, your usual command)
samtools fastq -T '*' monomers.bam \
  | minimap2 -ay -x map-hifi ref.fa - \
  | samtools view -b -o aligned.ns.bam

# 3. undo the cuts the enzyme never made, and get contacts
pore-c-aqb merge aligned.ns.bam \
    --output fragments.tsv.gz \
    --pairs contacts.pairs --sizes hg38.sizes.genome \
    --stats merge_stats.json --min-fragments 2

# 4. (diagnostic) which enzymes actually cut?
pore-c-aqb junctions aligned.ns.bam ref.fa --enzymes DpnII,NlaIII,HindIII
```

---

## Knowing where it is

A full library takes tens of minutes to digest and a 24 GB BAM a couple to
merge. Silence for that long is indistinguishable from a hang, so every long
command draws a bar on stderr:

```
digesting  [#########...............]  38%  1,504,220 reads  1,240/s  elapsed 12:07  eta 19:44
merging    [################........]  67%  139,204 concatemers  7,130/s  elapsed 0:19  eta 0:09
```

The percentage is **real**, not guessed. A BAM is BGZF-compressed and
`AlignmentFile.tell()` gives the position in the compressed file, so the
fraction of bytes consumed is a genuine measure of progress, and the ETA
follows from it. Reading from a pipe has no size to measure against, so the bar
falls back to a count and a rate rather than inventing a total.

It draws **only when stderr is a terminal**, so redirecting to a file or a
Nextflow log stays clean. `--progress` forces it on, `--no-progress` off.
`--quiet` does *not* silence it: that flag controls logging, and a bar is not a
log line.

---

## Fix 1 — several enzymes in one digest

A concatemer is digested **once**, by every enzyme in the tube at the same
time. Its fragments are delimited by the **union** of all recognition sites.

If the workflow is told about only one enzyme, a monomer that internally spans
a junction made by the *other* enzyme is never split. It contains two distinct
genomic loci, aligns to only one of them, and the rest is soft-clipped away.
**The contact is lost, silently, with no warning anywhere.**

```
 read:   ---------A---------+---------B---------+---------C---------
                          GATC                CATG
                         (DpnII)             (NlaIII)

 told only DpnII:   [----A----][-------- B + C --------]
                                         ^
                               two loci in one monomer; B aligns,
                               C is soft-clipped and thrown away

 told both:         [----A----][----B----][----C----]     all three kept
```

Measured end to end on real reads, through `samtools fastq` and `minimap2`:

| `--cutter` | monomers | aligned | soft-clipped bases |
|---|---:|---:|---:|
| `DpnII` | 17,076 | 7,390 | 218,363 |
| `DpnII,NlaIII` | 39,910 | 8,543 | **8,085** |

Soft-clipping falls by 96%. Those clipped bases *were* the second locus.

### Which enzymes can I use?

Any of the 729 Biopython knows a cut position for. `pore-c-aqb enzymes` shows
the ones that actually turn up in 3C protocols, and searches the rest by name
or by recognition site:

```
$ pore-c-aqb enzymes
  enzyme   site      cut (both strands)  sticky end  1 site per  recognition
  DpnII    GATC      N^GATC_N            5' GATC     256 bp      palindromic
  MboI     GATC      N^GATC_N            5' GATC     256 bp      palindromic
  NlaIII   CATG      _CATG^              3' CATG     256 bp      palindromic
  HinfI    GANTC     G^ANT_C             5' ANT      256 bp      palindromic
  ...
  NotI     GCGGCCGC  GC^GGCC_GC          5' GGCC     65,536 bp   palindromic

$ pore-c-aqb enzymes CATG      # by site
$ pore-c-aqb enzymes Dpn       # by name
$ pore-c-aqb enzymes --all     # all 729
```

The 334 enzymes with no defined cut position and the 25 that cut twice are left
out, because the digest rejects them anyway.

### How the union is computed

Each enzyme searches the read **separately**, then all positions go into one
set, are de-duplicated and sorted, and the read is cut once.

Searching once with a combined pattern is wrong, and measurably so: a plain
alternation (`GATC|CATG`) consumes its match and resumes after it, so
overlapping sites are lost. In `CATGATC` the `CATG` at 0 swallows the `GATC`
that starts at 3 — **243 wrong answers out of 1,600** when prototyped. Each
enzyme therefore gets its own pass, using a zero-width lookahead `(?=(GATC))`
so that even self-overlapping sites (`ATAT` in `ATATAT`) are all found.

A position cut by two enzymes is **one** cut, so the union is a set; the number
of such coincidences is reported separately, which is why the per-enzyme
percentages sum to slightly more than 100.

Adding an enzyme can only ever add cuts, never remove one — checked as a
property test, and on 16,758 real concatemers: **0 lost boundaries, 0 tiling
gaps**, 259,214 monomers becoming 634,479.

---

## Fix 2 — undo the cuts the enzyme never made

The digest cuts at every recognition site. The enzyme in the tube did not: real
digestion is incomplete, so a genuine restriction fragment usually contains
several uncut sites and ends up as several monomers that align head-to-tail.

Nothing warns you, and the damage is severe:

```
 reality:   locus X =====================  ligated to  ====== locus Y
                     (one uncut fragment)

 in silico: [x1][x2][x3][x4]                            [y1][y2]
                     |                                      |
                     +---------- 4 x 2 = 8 pairs -----------+
                          where there was ONE contact

            ...and x1-x2, x1-x3, x1-x4, x2-x3 ... all land on the DIAGONAL
```

On the library this tool came from: **11-fold duplication at 5 kb resolution,
18-fold at 25 kb, 44-76% of all pairs on the diagonal.** Juicer's
normalisations cannot repair it — only 14% of the bias is separable into a row
factor times a column factor, so a residual factor of 2.45 survives whatever
you normalise with. The distortion follows enzyme-site density, so it *deforms*
the map rather than merely scaling it.

### Why this cannot be fixed during the digest

From the read sequence alone, an uncut `GATC` inside a fragment and a `GATC`
reconstituted by ligation are **the same four letters**. There is no
information to tell them apart.

The information only appears after alignment: if the two pieces either side of
a site land next to each other on the genome, the site was never cut. That is
why this is a separate step, run on aligned monomers.

```bash
pore-c-aqb merge aligned.ns.bam --output fragments.tsv.gz \
    --pairs contacts.pairs --sizes hg38.sizes.genome --min-fragments 2
```

Two monomers merge when they are consecutive along the read **and** contiguous
on the genome: same chromosome, same strand, gap <= `--merge-gap` (100 bp by
default). It prints what it changed on *your* data:

```
Read 208,121 concatemers, 2,278,102 aligned monomers.
Merged into 450,213 fragments (5.06x fewer), keeping 129,864 molecules.

Cis contacts by separation, before and after merging:
  separation     before    after
  <1kb            25.1%     3.3%
  1-10kb          54.7%    23.1%
  10kb-1Mb        13.9%    47.8%
  >1Mb             6.2%    25.8%
```

Five cuts in six were false. The 10 kb - 1 Mb window, where TADs and loops
live, goes from 13.9% of pairs to 47.8%. The signal was there all along.

**A sanity check worth running on your own output.** Merged fragment length
should grow linearly with the number of monomers glued, because each uncut site
adds one restriction fragment. It does:

| monomers merged | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| median length | 356 bp | 747 bp | 1,170 bp | 1,592 bp | 1,998 bp |

**Order matters, and getting it wrong is silent.** Merging runs on the complete
chain of monomers; `--mapq` is applied *afterwards*, to the merged blocks.
Filter first and you break the chain — drop a middle monomer and its two
neighbours stop being adjacent, so they are counted as two loci in contact, a
contact that does not exist. During development this made the molecule count
*rise* with a stricter threshold (61,677 to 70,735), which is impossible and
was the giveaway. Two tests pin the ordering down.

Outputs: `fragments.tsv.gz` (one row per merged fragment, with how many
monomers were glued) and `contacts.pairs` (4DN v1.0, ready for `cooler cload
pairs` or `juicer_tools pre`).

---

## Fix 3 — which enzymes actually cut?

The digest report counts recognition **sites found in the reads**. That number
cannot tell a working enzyme from one that never left the freezer, because
every motif occurs in genomic DNA by chance — `AAGCTT` turns up every ~4 kb
whatever was in the tube.

> An earlier version of this tool labelled that column "cuts" and printed
> `HindIII 1,497 cuts (9.2%)` for an enzyme that demonstrably never cut. A user
> would have concluded the opposite of the truth. The column now says **sites
> found**, prints the chance rate beside the observed one, and states plainly
> that it is not evidence of cutting.

The test that *does* work needs aligned data. Inside a concatemer, the boundary
between two monomers that land far apart is a genuine ligation junction. If an
enzyme cut, its motif sits at those boundaries far more often than at random
positions on the same chromosomes:

```bash
pore-c-aqb junctions aligned.ns.bam ref.fa --enzymes DpnII,NlaIII,HindIII,HinfI
```

```
30,000 junction boundaries (MAPQ >= 20, jump > 1,000 bp)

  enzyme     site      at junctions   at random  enrichment   verdict
  DpnII      GATC             39.5%        1.3%       29.5x   yes, main enzyme
  NlaIII     CATG             10.3%        2.3%        4.4x   yes, secondary
  HindIII    AAGCTT            0.0%        0.1%        0.1x   no
  HinfI      GANTC             1.3%        1.5%        0.9x   no
```

That is the run this tool was written for: the protocol notes said HindIII,
`AAGCTT` sits *below* background at junctions, and the real second enzyme was
NlaIII. This command changes nothing on disk — it only reads.

Read the **enrichment** column, not the percentages: aligners soft-clip a few
bases at monomer ends, so the absolute rates understate. `--tol 10` recovers
DpnII to ~79% but inflates the random background too.

---

## Using it inside wf-pore-c

Three routes, from least to most invasive, in
[`docs/INTEGRATION.md`](docs/INTEGRATION.md). The simplest is to run the digest
yourself and feed the monomers to the rest of your pipeline.

To patch the workflow, [`patches/wf-pore-c-multicutter.patch`](patches/) applies
to a checkout of `epi2me-labs/wf-pore-c` and swaps the digest at **both** call
sites — `digest_align_annotate` invokes it once per chunking branch, with the
positional arguments in a different order each time:

```bash
pore-c-py digest "${meta.cutter}" ...                     # chunked: stdin
pore-c-py digest "concatemers.bam" "${meta.cutter}" ...   # not chunked
```

`pore-c-aqb` accepts **either order**, and implements the `--max_monomers`,
`--excluded_list`, `--excluded_bam`, `--recursive` and `--glob` options the
workflow passes. Patching one branch, or dropping those options, leaves the
workflow broken on one of its two paths.

---

## Correctness

Everything below is a test, not a claim.

**1. With one enzyme, output is identical to upstream.**

`tests/test_equivalence.py` runs the real `pore-c-py` inside the published
**wf-pore-c** image — the version the workflow actually uses, pore-c-py
**2.0.6** with biopython 1.82 — and diffs the BAMs. The strict comparison takes
every field and **every tag, including each tag's value type**, so a divergence
in something nobody thought to list still fails.

Validated on **16,758 real PacBio 3C reads, 259,214 monomers**, of which
203,342 carry base-modification tags: **zero differences**. The MM/ML
recomputation is vendored verbatim from upstream (`src/pore_c_aqb/_vendored.py`)
precisely so this holds.

> Pin the image, not the package name. The standalone `ontresearch/pore-c-py`
> image is a later 2.1.x which changed the mod-base tags (`ML` became a uint8
> array, `MN` was added). Diffing against it reports differences that are
> upstream's own version bump, not defects here.

**2. With several enzymes, the digest is a strict refinement.**

Every boundary of the single-enzyme digest survives, monomers tile each read
exactly, and reads reconstruct byte-for-byte from their monomers. Checked on
16,758 real concatemers and as a property test over random sequences.

**3. Merging is verified in both directions.**

Contiguous pieces merge; different chromosomes, different strands and distant
loci never do; reverse-strand geometry is handled; and two tests fail if the
MAPQ filter is moved before the merge.

**4. Awkward input fails with a sentence, not a traceback.**

`tests/test_robustness.py` feeds the commands empty BAMs, non-BAM files,
missing paths, unmapped-only reads, lowercase and `N`-containing sequences,
1-base reads, BAMs from other tools with no `Xc` tag, coordinate-sorted BAMs,
and extreme parameter values. Exit codes follow the convention: 0 success,
1 runtime error, 2 usage error.

**5. The progress bar never corrupts anything.**

`tests/test_progress.py` checks it is silent off a terminal, wiped on close and
on exception, redrawn at most a few times per burst, truncated rather than
wrapped, and that it reports no percentage when it cannot know one.

```
263 tests | 96% coverage | biopython 1.82 and 1.88 | flake8 clean
```

### A reproducibility bug found along the way

`Bio.Restriction.search()` changed its behaviour at sequence boundaries:

```python
NlaIII.search(Seq("CATGCCTTCTGTGCGAGCCC"))
# biopython 1.82 -> [5]     (the version inside wf-pore-c)
# biopython 1.88 -> []
```

1.88 drops sites whose second-strand cut lands exactly on a sequence edge —
always at position 0 or the very end, always producing a few-base monomer that
never aligns. Rare (2 cases in 350 pairs built to hit boundaries, 0 in 2,400
random sequences), but **deterministic**: the same BAM gives different monomer
counts on two machines, which a tool meant to be cited cannot do.

`src/pore_c_aqb/sites.py` therefore locates sites itself for palindromic
enzymes — every enzyme used in 3C/Hi-C — under one documented rule. Verified
over **2,400 comparisons across 20 palindromic enzymes: zero disagreements**,
and 634,479 identical monomers from the same real library under both versions.

Non-palindromic Type IIS enzymes (BsrI, AcuI, ...) cut at a *distance* from
their site; the reverse-strand geometry is fiddly, an early attempt was off by
2 bases on every one of them, and none is used to build a Hi-C library. They
are delegated to Biopython rather than shipped on under-tested arithmetic. The
principle throughout: *use our code where it has been proven, delegate where it
has not.*

Enzymes Biopython cannot digest safely are **rejected up front** with a clear
message: 334 with no defined cut position, 25 that cut twice.

---

## Performance

16,758 real 3C reads (mean 4.4 kb), single core, second run of two:

| | time | throughput |
|---|---:|---:|
| `pore-c-py digest DpnII` (upstream, in its container) | 42.3 s | 395 reads/s |
| `pore-c-aqb digest DpnII` | 39.9 s | 420 reads/s |
| `pore-c-aqb digest DpnII,NlaIII` | 60.8 s | 275 reads/s |

No regression against upstream. A second enzyme costs about 50% more, split
between the extra site search and the 2.4x larger output. The digest is not the
bottleneck — alignment of the same reads takes an order of magnitude longer.

---

## Documentation

- [`docs/DESIGN.md`](docs/DESIGN.md) — every design choice, with the
  measurement behind it, including the ideas that were tried and rejected
- [`docs/INTEGRATION.md`](docs/INTEGRATION.md) — wiring it into wf-pore-c
- [`CHANGELOG.md`](CHANGELOG.md)
- [`NOTICE`](NOTICE) — exactly which files came from upstream, and every change

---

## Licence and attribution

Derived from [`pore-c-py`](https://github.com/epi2me-labs/pore-c-py) 2.0.6 by
Oxford Nanopore Technologies PLC, and distributed under the same **Oxford
Nanopore Technologies PLC Public License Version 1.0** (see [`LICENSE`](LICENSE)).

That licence restricts use to **Research Purposes**. It is not an OSI-approved
open-source licence; read it before depending on this in a commercial setting.

This project is not affiliated with or endorsed by Oxford Nanopore
Technologies.
