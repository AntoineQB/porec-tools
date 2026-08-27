# Design notes

Every choice below was made against a measurement. The numbers are reproducible
with the commands given.

---

## 1. Union of cut points, not one digest per enzyme

A concatemer is digested **once**, by all enzymes present simultaneously. Its
fragments are delimited by the union of every recognition site.

Running the digest once per enzyme and merging the outputs afterwards would be
wrong in kind, not just in degree: it produces *overlapping monomer sets*
rather than one partition of the read. A read cut by DpnII into `[0,50)[50,90)`
and by NlaIII into `[0,30)[30,90)` has a true double digest of
`[0,30)[30,50)[50,90)`, which neither single digest contains.

So `find_cut_points` takes the union of positions, then reuses upstream's
`splits_to_intervals` unchanged.

**De-duplication.** Two enzymes can cut at the same base. That is one cut in
the tube, so the union is a `set`. The count of such coincidences is reported
separately (`shared_cut_points`) so the per-enzyme percentages still add up
honestly.

---

## 2. Biopython for searching, not a hand-rolled scanner, then partly the
   opposite

### First attempt: one regex pass over all motifs

The obvious optimisation is a single alternation regex
(`GATC|CATG|GA[ACGT]TC|...`) scanned once, instead of one Biopython pass per
enzyme. Prototyped and measured:

| enzymes | regex | Biopython | speed-up |
|---|---:|---:|---:|
| 1 | 0.028 ms | 0.105 ms | 3.7x |
| 2 | 0.180 ms | 0.223 ms | 1.2x |
| 4 | 0.289 ms | 0.414 ms | 1.4x |

**And 243 wrong answers out of 1,600.** A plain alternation scan consumes its
match and resumes after it, so overlapping sites are lost: in `CATGATC` the
`CATG` match at 0 swallows the `GATC` that starts at 3.

For 1.2 to 1.4x on the realistic case, against a demonstrated correctness bug, the
trade was refused. The digest is not the bottleneck anyway. Alignment costs an
order of magnitude more.

### Then: a reproducibility problem forced a scanner after all

```python
NlaIII.search(Seq("CATGCCTTCTGTGCGAGCCC"))
# biopython 1.82 -> [5]
# biopython 1.88 -> []
```

Biopython 1.88 drops sites whose second-strand cut lands exactly on a sequence
edge; 1.82 keeps them. Measured over 350 enzyme/sequence pairs: **2
disagreements, both at position 0 or at the very end**.

The monomers involved are a few bases long and never align, so no result
changes. But the *same BAM digested on two machines yields different monomer
counts*, and a tool meant to be cited cannot do that.

`sites.py` therefore scans for sites itself, using a **lookahead** regex,
`(?=(GATC))`, which is zero-width and so reports overlapping matches, avoiding
the bug that killed the first attempt. One documented rule: a cut at
`site_start + fst5` for every occurrence, kept when inside `[0, len]`.

Validation: **2,400 comparisons across 20 palindromic enzymes, zero
disagreements outside the documented boundary case.**

### Except for Type IIS enzymes

Non-palindromic enzymes cut at a *distance* from their site:

```
BsrI   ACTG_GN^N                          fst5=6  for a 5 bp site
AcuI   CTGAAGNNNNNNNNNNNNNN_NN^N          fst5=22 for a 6 bp site
```

Getting the reverse-strand geometry right for these is fiddly, an early attempt
was off by 2 bases on every one of them, and **none is used to build a Hi-C
library**. They are delegated to Biopython rather than shipped on under-tested
arithmetic. The version-dependent edge case can still bite there; that is
stated in the module docstring rather than hidden.

This is the general principle applied throughout: *use our code where it has
been proven, delegate where it has not.*

---

## 3. Vendoring the MM/ML recomputation verbatim

`get_subread_modified_bases` rebuilds the methylation tags for each monomer
after trimming. It is the subtlest code in the digest: an off-by-one silently
corrupts every downstream methylation analysis, and nothing would fail loudly.

It is copied **unchanged** into `_vendored.py`, with attribution. Any
"improvement" would break the byte-identity guarantee, which
`tests/test_equivalence.py` enforces by diffing against the real upstream tool.

Validated against the real tool in the wf-pore-c container (pore-c-py 2.0.6)
on 16,758 real PacBio reads: 259,214 monomers, 203,342 of them carrying
base-modification tags, **zero differences on every field and every tag,
value types included**.

Two traps here, both hit during development:

* Compare against the **right version**. The standalone `ontresearch/pore-c-py`
  image is a later 2.1.x that writes `ML` as a uint8 array and adds `MN`;
  2.0.6 writes `ML` as a string and has no `MN`. Diffing against 2.1.x
  produced three "defects" that were nothing of the sort. wf-pore-c ships
  2.0.6, so 2.0.6 is what this tracks.
* Compare **everything**. The first version of the comparison listed the tags
  it cared about (`MI`, `Xc`, `MM`, `ML`) and ran on synthetic reads carrying
  none of them, so it would have passed with the mod-base handling entirely
  broken. It now compares every tag *and its value type*, on reads that
  actually carry `MM`/`ML`. Deliberately writing `ML` as a string makes three
  tests fail, which is the property a test of this kind needs.

---

## 4. Validating enzymes before reading any data

A typo in the second of three enzymes used to surface after the first read was
processed. `resolve_enzymes` resolves and checks every name first, so a bad
specification fails in milliseconds rather than after hours.

Also handled at resolution time:

- **Isoschizomers collapsed.** `DpnII,MboI` both recognise `GATC` and cut
  identically; searching twice is wasted work and would double-count in the
  per-enzyme report.
- **Rejected**: enzymes that cut twice (two cut points per site, which
  `splits_to_intervals` does not model), and enzymes with no defined cut
  position, for instance `HpyUM037X`, whose site is `TNGGNAG|GTGGNAG`.
- **Typos** get close-match suggestions, case-insensitive first, so `dpnii`
  points at `DpnII`.

---

## 5. Per-enzyme statistics, and what they are not

Not decoration. The dataset that motivated this tool had been digested with
DpnII alone, while a fifth of its ligation junctions sat on `CATG`. Finding
that out required a bespoke junction analysis after the fact.

`DigestStats` makes the per-enzyme contribution a line of routine output:

```
  enzyme   site    sites found  % of cuts  1 per (observed)  1 per (chance)
  DpnII    GATC         14,721      37.5%            357 bp          256 bp
  NlaIII   CATG         23,104      58.8%            227 bp          256 bp
  HindIII  AAGCTT        1,497       3.8%          3,507 bp        4,096 bp
```

### The mistake this section exists to record

A pre-release version labelled that column **cuts**, and the caveat in the
module docstring called it a way to answer "did this enzyme actually cut?".
Both were wrong, and dangerously so. On the very run above it printed:

```
  HindIII      AAGCTT            1,497 cuts (  9.2%)
```

HindIII did not cut this library. Junction analysis puts `AAGCTT` *below*
background at ligation junctions. What the number counts is occurrences of the
motif in the reads, and `AAGCTT` occurs every ~4 kb in human DNA no matter what
was in the tube. A user reading "9.2% cuts" would conclude the enzyme worked.
That is the opposite of the truth, produced by the tool's own summary line.

Two things were changed:

1. **The column says `sites found`,** and the chance rate is printed beside the
   observed one so the reader can see they match. A closing note states that
   the number is not proof of cutting.
2. **The real test ships as a second command,** `pore-c-aqb-junctions`, which
   measures motif enrichment at ligation junctions against a random background
   on the same chromosomes. It needs alignments, which is exactly why the
   digest cannot answer the question.

What the digest table legitimately gives you: confirmation the enzyme was
applied, a typo check (`0 sites` warns), and the share of fragmentation each
enzyme contributes, the number you need to decide whether a second enzyme
earned its place.

### A diagnostic that was tried and rejected

Before settling on this, the read termini were tried: a concatemer's outer ends
are genuine cuts, so a read starting on an enzyme's cut should begin with
`site[fst5:]`. On 1,036 real reads it returned **0.0% for every enzyme,
including DpnII**, which certainly cut. PacBio reads carry adapters and
barcodes at their ends, so the termini are not cut sites. Recording it here so
it is not re-attempted.

---

## 6. Performance

16,758 real 3C reads, mean 4.4 kb, single core, second run of two:

| | time | throughput |
|---|---:|---:|
| upstream `pore-c-py digest DpnII` | 42.3 s | 395 reads/s |
| `pore-c-aqb digest DpnII` | 39.9 s | 420 reads/s |
| `pore-c-aqb digest DpnII,NlaIII` | 60.8 s | 275 reads/s |

No regression. The second enzyme costs ~50%, split between the extra search and
a 2.4x larger output to write (259,214 -> 634,479 monomers).

`--threads` only affects BAM compression (14.1 s vs 15.0 s at 4 threads); the
work is dominated by MM/ML recomputation, which is per-monomer and vendored.

Extrapolated: a 1M-read library takes ~1 h with two enzymes. Acceptable for a
step run once per dataset, and small beside alignment.

### Reproducing these numbers

```bash
pytest                                    # correctness
python -m pore_c_aqb.cli digest DpnII,NlaIII reads.bam --output /dev/null --stats s.tsv
```

The equivalence tests against upstream need Docker and the
`ontresearch/wf-pore-c` image; they skip cleanly without it.
