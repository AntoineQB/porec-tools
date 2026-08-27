# porec-tools

**A modified version of Oxford Nanopore's [`pore-c-py`](https://github.com/epi2me-labs/pore-c-py) / [`wf-pore-c`](https://github.com/epi2me-labs/wf-pore-c), for Pore-C and PacBio CiFi libraries.**

Not a new pipeline: this is the upstream digest, kept byte-for-byte identical
where it was already right, with three things fixed that were silently costing
contacts:

| | Problem | Fix |
|---|---|---|
| 1 | `--cutter` accepts one enzyme. Many protocols use two or more, and the second enzyme's junctions are never cut. | `porec digest DpnII,NlaIII` |
| 2 | The *in silico* digest cuts at every site, including the many the enzyme never cut. This inflates the Hi-C diagonal with contacts that never happened. | `porec merge` |
| 3 | Nothing tells you which enzymes actually cut your library. | `porec junctions` |

Plus one reproducibility bug in a dependency, described [below](#a-reproducibility-bug-found-along-the-way).

With a single enzyme, the digest output is identical to upstream. It is verified
record by record, tag by tag, against the real tool. See [Correctness](#correctness).

<details>
<summary><b>Everything that changed, at a glance</b></summary>

| Change | Where | Runs on |
|---|---|---|
| `--cutter` takes a list of enzymes; cut points are their union | `digest` | unaligned |
| Enzyme names validated up front, with typo suggestions | `digest` | either |
| Enzymes with no defined cut (334) or that cut twice (25) rejected, not mis-digested | `digest` | either |
| Site search made independent of the Biopython version | `digest` | either |
| `--stats` : per-enzyme site report, as TSV | `digest` | either |
| `--dry-run` : show what each enzyme will do, reading no data | `digest` | either |
| `@PG` provenance line in the output BAM header | `digest` | either |
| Both positional orders accepted (`ENZYME INPUT` and `INPUT ENZYME`) | `digest` | either |
| `--max_monomers`, `--excluded_list`, `--excluded_bam`, `--recursive`, `--glob` implemented, as the workflow passes them | `digest` | either |
| **New command:** merge fragments the enzyme never cut, write 4DN `.pairs` | `merge` | **aligned** |
| **New command:** motif enrichment at ligation junctions | `junctions` | **aligned** |
| **New command:** list and search usable enzymes | `enzymes` | either |
| Progress bar with a real percentage and ETA | all long commands | either |
| Clear errors instead of tracebacks on bad input | all | either |

Full detail in [`NOTICE`](NOTICE), which lists every file taken from upstream
and every change made to it.

</details>

---

## Where it fits

Two of the three fixes act on **unaligned** reads, one on **aligned** ones.
That split is not a design preference. It is forced by what information exists
at each point, and it is the thing to understand before using the tool.

```
   concatemers.bam   (unaligned reads, straight off the sequencer)
          |
          v
   +---------------------------+
   |  porec digest        |  cut each read into monomers
   |  DpnII,NlaIII             |  <- FIX 1
   +---------------------------+
          |  monomers, still unaligned
          v
   +---------------------------+
   |  minimap2                 |  NOW each monomer gets a genomic position.
   |  (unchanged)              |  Only here does it become knowable which
   +---------------------------+  cuts were real.
          |  aligned.ns.bam   (name-sorted: monomers of a read stay together)
          |
          +--------------------------------+
          v                                v
   +---------------------------+   +---------------------------+
   |  porec merge         |   |  porec junctions     |
   |  undo the false cuts      |   |  which enzymes cut?       |
   |  <- FIX 2                 |   |  <- FIX 3 (read-only)     |
   +---------------------------+   +---------------------------+
          |  fragments.tsv.gz + contacts.pairs
          v
      cooler / juicer
```

**Why the line falls there.** In an unaligned read, a `GATC` that the enzyme
never cut and a `GATC` rebuilt by ligation are the same four letters. Nothing
distinguishes them. Only once each monomer has a genomic position can you see
that two of them land side by side, which means the site between them was
never cut. So `merge` and `junctions` *cannot* run before alignment, and
`digest` cannot wait for it.

### Where does alignment happen in wf-pore-c?

Inside the same process as the digest. The workflow streams the whole thing
through one pipe, with no intermediate file:

```bash
pore-c-py digest "concatemers.bam" "${meta.cutter}" ... |
samtools fastq -T '*' |
minimap2 -ay -t N reference.fasta.mmi - |
pore-c-py annotate - "${meta.alias}" --monomers --stdout |
tee "${meta.alias}_out.ns.bam" |
samtools sort -o "${meta.alias}.cs.bam" -
```

So `porec merge` does not go inside that pipe. You run it afterwards, on
the `*_out.ns.bam` the workflow already writes:

```
   wf-pore-c  (digest | fastq | minimap2 | annotate)  ->  SAMPLE_out.ns.bam
                                                                  |
                                          porec merge  <-----+
```

Use the **name-sorted** `.ns.bam`, not the coordinate-sorted `.cs.bam`: merging
needs the monomers of one read to be adjacent. The tool refuses a
coordinate-sorted file rather than reading it wrongly.

---

## Install

```bash
pip install git+https://github.com/YOUR-USERNAME/porec-tools.git
```

> **Before the first push**, point the URLs at your own account:
> ```bash
> grep -rl YOUR-USERNAME . --exclude-dir=.git \
>   | xargs sed -i 's/YOUR-USERNAME/your-github-handle/g'
> ```

Needs Python >= 3.8, `pysam` and `biopython`. The last two are pulled in
automatically.
Three commands are installed:

```bash
porec enzymes      # which enzymes can I use?
porec digest
porec merge        # also available as porec-merge
porec junctions    # also available as porec-junctions
```

For development:

```bash
git clone https://github.com/YOUR-USERNAME/porec-tools.git
cd porec-tools
pip install -e ".[dev]"
pytest
```

---

## Quick start

```bash
# 0. find your enzymes, then check what the digest will do with them
porec enzymes                    # the 3C/Hi-C shortlist
porec enzymes CATG               # everything cutting CATG
porec digest DpnII,NlaIII --dry-run

# 1. digest with both enzymes
porec digest DpnII,NlaIII reads.bam \
    --output monomers.bam --stats digest_stats.tsv --threads 8

# 2. align  <-- NOTHING OF THIS TOOL RUNS HERE. Your usual command, unchanged.
#              Everything below needs the genomic positions it produces.
samtools fastq -T '*' monomers.bam \
  | minimap2 -ay -x map-hifi ref.fa - \
  | samtools view -b -o aligned.ns.bam        # name-sorted, NOT coordinate

# 3. undo the cuts the enzyme never made, and get contacts  (needs step 2)
porec merge aligned.ns.bam \
    --output fragments.tsv.gz \
    --pairs contacts.pairs --sizes hg38.sizes.genome \
    --stats merge_stats.json --min-fragments 2

# 4. (diagnostic, needs step 2) which enzymes actually cut?
porec junctions aligned.ns.bam ref.fa --enzymes DpnII,NlaIII,HindIII
```

---

## Progress: knowing how far along a run is

A full library takes tens of minutes to digest and a 24 GB BAM a couple to
merge. With no output during that time you cannot tell a slow run from a
stuck one, so the long commands draw a bar on stderr:

```
digesting  [#########...............]  38%  1,504,220 reads  1,240/s  elapsed 12:07  eta 19:44
merging    [################........]  67%  139,204 concatemers  7,130/s  elapsed 0:19  eta 0:09
```

The percentage is real, not guessed. A BAM is BGZF-compressed and
`AlignmentFile.tell()` gives the position in the compressed file, so the
fraction of bytes consumed is a genuine measure of progress, and the ETA
follows from it. Reading from a pipe has no size to measure against, so the bar
falls back to a count and a rate rather than inventing a total.

It draws only when stderr is a terminal, so redirecting to a file or a
Nextflow log stays clean. `--progress` forces it on, `--no-progress` off.
`--quiet` does *not* silence it: that flag controls logging, and a bar is not a
log line.

---

## Fix 1. Several enzymes in one digest

A concatemer is digested once, by every enzyme in the tube at the same
time. Its fragments are delimited by the union of all recognition sites.

If the workflow is told about only one enzyme, a monomer that internally spans
a junction made by the *other* enzyme is never split. It contains two distinct
genomic loci, aligns to only one of them, and the rest is soft-clipped away.
The contact is lost, with no warning anywhere.

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
| `DpnII,NlaIII` | 39,910 | 8,543 | 8,085 |

Soft-clipping falls by 96%. Those clipped bases *were* the second locus.

### Which enzymes can I use?

Any of the 729 Biopython knows a cut position for. `porec enzymes` shows
the ones that actually turn up in 3C protocols, and searches the rest by name
or by recognition site:

```
$ porec enzymes
  enzyme   site      cut (both strands)  sticky end  1 site per  recognition
  DpnII    GATC      N^GATC_N            5' GATC     256 bp      palindromic
  MboI     GATC      N^GATC_N            5' GATC     256 bp      palindromic
  NlaIII   CATG      _CATG^              3' CATG     256 bp      palindromic
  HinfI    GANTC     G^ANT_C             5' ANT      256 bp      palindromic
  ...
  NotI     GCGGCCGC  GC^GGCC_GC          5' GGCC     65,536 bp   palindromic

$ porec enzymes CATG      # by site
$ porec enzymes Dpn       # by name
$ porec enzymes --all     # all 729
```

The 334 enzymes with no defined cut position and the 25 that cut twice are left
out, because the digest rejects them anyway.

### How the union is computed

Each enzyme searches the read separately, then all positions go into one
set, are de-duplicated and sorted, and the read is cut once.

Searching once with a combined pattern is wrong, and measurably so: a plain
alternation (`GATC|CATG`) consumes its match and resumes after it, so
overlapping sites are lost. In `CATGATC` the `CATG` at 0 swallows the `GATC`
that starts at 3. When prototyped this gave **243 wrong answers out of
1,600**. Each
enzyme therefore gets its own pass, using a zero-width lookahead `(?=(GATC))`
so that even self-overlapping sites (`ATAT` in `ATATAT`) are all found.

A position cut by two enzymes is one cut, so the union is a set; the number
of such coincidences is reported separately, which is why the per-enzyme
percentages sum to slightly more than 100.

Adding an enzyme can only ever add cuts, never remove one. That is checked as
a property test, and on 16,758 real concatemers: **0 lost boundaries, 0 tiling
gaps**, 259,214 monomers becoming 634,479.

---

## Fix 2. Undo the cuts the enzyme never made

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
normalisations cannot repair it, because only 14% of the bias is separable into a row
factor times a column factor, so a residual factor of 2.45 survives whatever
you normalise with. The distortion follows enzyme-site density, so it *deforms*
the map rather than just scaling it.

### Why this cannot be fixed during the digest

From the read sequence alone, an uncut `GATC` inside a fragment and a `GATC`
reconstituted by ligation are the same four letters. There is no
information to tell them apart.

The information only appears after alignment: if the two pieces either side of
a site land next to each other on the genome, the site was never cut. That is
why this is a separate step, run on aligned monomers.

```bash
porec merge aligned.ns.bam --output fragments.tsv.gz \
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
live, goes from 13.9% of pairs to 47.8%.

A sanity check you can run on your own output: Merged fragment length
should grow linearly with the number of monomers glued, because each uncut site
adds one restriction fragment. It does:

| monomers merged | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| median length | 356 bp | 747 bp | 1,170 bp | 1,592 bp | 1,998 bp |

Order matters here, and getting it wrong is silent. Merging runs on the complete
chain of monomers; `--mapq` is applied *afterwards*, to the merged blocks.
Filter first and you break the chain. Drop a middle monomer and its two
neighbours stop being adjacent, so they are counted as two loci in contact, a
contact that does not exist. During development this made the molecule count
*rise* with a stricter threshold (61,677 to 70,735), which is impossible and
is how the ordering bug was found. Two tests cover it.

Outputs: `fragments.tsv.gz` (one row per merged fragment, with how many
monomers were glued) and `contacts.pairs` (4DN v1.0, ready for `cooler cload
pairs` or `juicer_tools pre`).

---

## Fix 3. Which enzymes actually cut?

The digest report counts recognition sites found in the reads. That number
cannot tell a working enzyme from one that never left the freezer, because
every motif occurs in genomic DNA by chance. `AAGCTT` turns up every ~4 kb
whatever was in the tube.

> An earlier version of this tool labelled that column "cuts" and printed
> `HindIII 1,497 cuts (9.2%)` for an enzyme that demonstrably never cut. A user
> would have concluded that it worked. The column now says **sites found**, prints the chance rate beside the observed one, and states plainly
> that it is not evidence of cutting.

The test that *does* work needs aligned data. Inside a concatemer, the boundary
between two monomers that land far apart is a genuine ligation junction. If an
enzyme cut, its motif sits at those boundaries far more often than at random
positions on the same chromosomes:

```bash
porec junctions aligned.ns.bam ref.fa --enzymes DpnII,NlaIII,HindIII,HinfI
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
NlaIII. This command changes nothing on disk, it only reads.

Read the enrichment column, not the percentages: aligners soft-clip a few
bases at monomer ends, so the absolute rates understate. `--tol 10` recovers
DpnII to ~79% but inflates the random background too.

---

## Command reference

Every option of every command. `--help` on any of them prints the same, with
more explanation.

### `porec enzymes`, what can I digest with?

```
porec enzymes [QUERY] [--all]

  QUERY          A name fragment ('Dpn') or a recognition site ('GATC').
                 Searching implies --all.
  --all          Every usable enzyme (729), not just the 3C/Hi-C shortlist.
```

*Entirely new, no upstream equivalent.*

### `porec digest`, cut reads into monomers (unaligned input)

Drop-in for `pore-c-py digest`: every upstream option is accepted with the same
meaning, and both positional orders work, because wf-pore-c uses both.

```
porec digest ENZYME [INPUT ...]   [options]
porec digest [INPUT ...] ENZYME   [options]      # upstream's order

  ENZYME               One or more enzyme names. 'DpnII', 'DpnII,NlaIII',
                       'DpnII+NlaIII' all work.
  INPUT                Unaligned BAM(s), a directory, or '-' for stdin
                       (the default).

  --output PATH        Output BAM ('-' for stdout).                    [-]
  --header PATH        BAM whose header to copy. Required with stdin.
  --stats PATH         Per-enzyme site report, as TSV.                  NEW
  --dry-run            Print what each enzyme will do, then exit.       NEW
                       Reads no data, so it costs nothing before a
                       long run.
  --remove_tags ...    Extra SAM tags to strip.
  --max_reads N        Take only the first N concatemers.
  --max_monomers N     Drop a concatemer cut into more than N monomers.
  --excluded_list PATH Names of the reads dropped by --max_monomers.
  --excluded_bam PATH  The dropped reads themselves.
  --recursive          Search an input directory recursively.
  --glob PATTERN       Which files to take from a directory.       [*.bam]
  --threads N          Threads for BAM compression.                    [1]
  --debug / --quiet    Logging verbosity.
  --logfile PATH       Also write logs here.
  --progress / --no-progress                                           NEW
  --version                                                            NEW
```

Also new, without a flag: several enzymes at once, names validated before a
single read is touched, and a `@PG` provenance line in the output header.

### `porec merge`, undo the false cuts (ALIGNED input)

**Runs after alignment.** See [Where it fits](#where-it-fits).

```
porec merge ALIGNED_BAM [options]

  ALIGNED_BAM          Aligned monomer BAM, grouped by read name. This is
                       the workflow's *.ns.bam. A coordinate-sorted BAM is
                       rejected rather than silently misread.

  --output PATH        Merged fragments, as TSV. '.gz' compresses.     [-]
  --merge-gap BP       Largest gap between two monomers still counted
                       as one uncut fragment.                        [100]
  --mapq Q             Minimum mapping quality, applied AFTER merging.  [1]
  --min-fragments N    Keep only molecules with at least N fragments.
                       Use 2 for contacts, 3 for multi-way analysis.   [1]
  --min-length BP      Drop merged fragments shorter than this.         [0]
  --pairs PATH         Also write contacts as 4DN .pairs.
  --sizes PATH         Chromosome sizes, to order the .pairs upper
                       triangle consistently with your genome.
  --min-sep BP         Drop cis pairs closer than this. Blunt; the
                       rigorous criterion is the number of restriction
                       sites between the two fragments.                 [0]
  --stats PATH         Summary as JSON, with every parameter recorded.
  --quiet              Print nothing but errors.
  --progress / --no-progress
  --version
```

*Entirely new, no upstream equivalent.*

### `porec junctions`, which enzymes actually cut? (ALIGNED input)

Runs after alignment, and writes nothing: it only reads.

```
porec junctions ALIGNED_BAM REFERENCE --enzymes LIST [options]

  ALIGNED_BAM          Name-sorted aligned monomer BAM.
  REFERENCE            Indexed reference FASTA (.fai required).
  --enzymes LIST       Comma-separated, e.g. 'DpnII,NlaIII,HindIII'.

  --mapq Q             Minimum monomer mapping quality.                [20]
  --min-jump BP        Minimum genomic distance for a boundary to count
                       as a ligation junction rather than an uncut
                       site.                                        [1000]
  --tol BP             Slack between the motif's cut and the alignment
                       boundary; covers the overhang and a little
                       soft-clipping.                                   [2]
  --max-junctions N    Stop after this many boundaries (0 = all).  [200000]
  --seed N             Seed for the random background, so a run is
                       reproducible.                                    [0]
  --progress / --no-progress
  --version
```

*Entirely new, no upstream equivalent.*

---

## Using it inside wf-pore-c

Three routes, from least to most invasive, in
[`docs/INTEGRATION.md`](docs/INTEGRATION.md). The simplest is to run the digest
yourself and feed the monomers to the rest of your pipeline.

To patch the workflow, [`patches/wf-pore-c-multicutter.patch`](patches/) applies
to a checkout of `epi2me-labs/wf-pore-c` and swaps the digest at both call
sites. `digest_align_annotate` invokes it once per chunking branch, with the
positional arguments in a different order each time:

```bash
pore-c-py digest "${meta.cutter}" ...                     # chunked: stdin
pore-c-py digest "concatemers.bam" "${meta.cutter}" ...   # not chunked
```

`porec` accepts either order, and implements the `--max_monomers`,
`--excluded_list`, `--excluded_bam`, `--recursive` and `--glob` options the
workflow passes. Patching one branch, or dropping those options, leaves the
workflow broken on one of its two paths.

---

## Correctness

Everything below is a test, not a claim.

**1. With one enzyme, output is identical to upstream.**

`tests/test_equivalence.py` runs the real `pore-c-py` inside the published
wf-pore-c image, which is the version the workflow actually uses (pore-c-py
**2.0.6** with biopython 1.82), and diffs the BAMs. The strict comparison takes
every field and every tag, including each tag's value type, so a divergence
in something nobody thought to list still fails.

Validated on **16,758 real PacBio 3C reads, 259,214 monomers**, of which
203,342 carry base-modification tags: zero differences. The MM/ML
recomputation is vendored verbatim from upstream (`src/porec_tools/_vendored.py`)
so that this holds.

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

1.88 drops sites whose second-strand cut lands exactly on a sequence edge,
always at position 0 or the very end, always producing a few-base monomer that
never aligns. Rare (2 cases in 350 pairs built to hit boundaries, 0 in 2,400
random sequences), but deterministic: the same BAM gives different monomer
counts on two machines, which a tool meant to be cited cannot do.

`src/porec_tools/sites.py` therefore locates sites itself for palindromic
enzymes, which covers every enzyme used in 3C/Hi-C, under one documented rule. Verified
over **2,400 comparisons across 20 palindromic enzymes: zero disagreements**,
and 634,479 identical monomers from the same real library under both versions.

Non-palindromic Type IIS enzymes (BsrI, AcuI, ...) cut at a *distance* from
their site; the reverse-strand geometry is fiddly, an early attempt was off by
2 bases on every one of them, and none is used to build a Hi-C library. They
are delegated to Biopython rather than shipped on under-tested arithmetic. The
same approach is used throughout: our own code where it has been checked
against a reference, Biopython where it has not.

Enzymes Biopython cannot digest safely are rejected up front with a clear
message: 334 with no defined cut position, 25 that cut twice.

---

## Performance

16,758 real 3C reads (mean 4.4 kb), single core, second run of two:

| | time | throughput |
|---|---:|---:|
| `pore-c-py digest DpnII` (upstream, in its container) | 42.3 s | 395 reads/s |
| `porec digest DpnII` | 39.9 s | 420 reads/s |
| `porec digest DpnII,NlaIII` | 60.8 s | 275 reads/s |

No regression against upstream. A second enzyme costs about 50% more, split
between the extra site search and the 2.4x larger output. The digest is not the
bottleneck. Alignment of the same reads takes an order of magnitude longer.

---

## Documentation

- [`docs/DESIGN.md`](docs/DESIGN.md), every design choice with the
  measurement behind it, including the ideas that were tried and rejected
- [`docs/INTEGRATION.md`](docs/INTEGRATION.md), wiring it into wf-pore-c
- [`CHANGELOG.md`](CHANGELOG.md)
- [`NOTICE`](NOTICE), exactly which files came from upstream and every change

---

## Licence and attribution

Derived from [`pore-c-py`](https://github.com/epi2me-labs/pore-c-py) 2.0.6 by
Oxford Nanopore Technologies PLC, and distributed under the same **Oxford
Nanopore Technologies PLC Public License Version 1.0** (see [`LICENSE`](LICENSE)).

That licence restricts use to **Research Purposes**. It is not an OSI-approved
open-source licence; read it before depending on this in a commercial setting.

This project is not affiliated with or endorsed by Oxford Nanopore
Technologies.
