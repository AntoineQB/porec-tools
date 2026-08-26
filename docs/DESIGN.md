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
`[0,30)[30,50)[50,90)` — which neither single digest contains.

So `find_cut_points` takes the union of positions, then reuses upstream's
`splits_to_intervals` unchanged.

**De-duplication.** Two enzymes can cut at the same base. That is one cut in
the tube, so the union is a `set`. The count of such coincidences is reported
separately (`shared_cut_points`) so the per-enzyme percentages still add up
honestly.

---

## 2. Biopython for searching, not a hand-rolled scanner — then partly the
   opposite

### First attempt: one regex pass over all motifs

The obvious optimisation is a single alternation regex
(`GATC|CATG|GA[ACGT]TC|…`) scanned once, instead of one Biopython pass per
enzyme. Prototyped and measured:

| enzymes | regex | Biopython | speed-up |
|---|---:|---:|---:|
| 1 | 0.028 ms | 0.105 ms | 3.7× |
| 2 | 0.180 ms | 0.223 ms | 1.2× |
| 4 | 0.289 ms | 0.414 ms | 1.4× |

**And 243 wrong answers out of 1,600.** A plain alternation scan consumes its
match and resumes after it, so overlapping sites are lost: in `CATGATC` the
`CATG` match at 0 swallows the `GATC` that starts at 3.

For 1.2–1.4× on the realistic case, against a demonstrated correctness bug, the
trade was refused. The digest is not the bottleneck anyway — alignment costs an
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

`sites.py` therefore scans for sites itself — using a **lookahead** regex,
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

Validated on 4,227 real PacBio reads: 65,456 monomers, **zero differing
records**, including 51,332 carrying `MM` tags.

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
  position — for instance `HpyUM037X`, whose site is `TNGGNAG|GTGGNAG`.
- **Typos** get close-match suggestions, case-insensitive first, so `dpnii`
  points at `DpnII`.

---

## 5. Per-enzyme statistics

Not decoration. The dataset that motivated this tool had been digested with
DpnII alone, while a fifth of its ligation junctions sat on `CATG`. Finding
that out required a bespoke junction analysis after the fact.

`DigestStats` makes it a line of routine output:

```
  DpnII        GATC             61,229 cuts ( 39.1%)
  NlaIII       CATG             95,424 cuts ( 61.0%)
  HindIII      AAGCTT                0 cuts (  0.0%)   <-- no sites found, check the protocol
```

An enzyme that failed in the reaction, or was named by mistake, is visible in
the first minute rather than three analysis steps later.

---

## 6. Performance

4,227 real 3C reads, mean 4.4 kb, single core:

| | time | throughput |
|---|---:|---:|
| upstream `pore-c-py digest DpnII` | ~10.8 s | 390 reads/s |
| `pore-c-aqb digest DpnII` | ~10.0 s | 430 reads/s |
| `pore-c-aqb digest DpnII,NlaIII` | ~15.0 s | 282 reads/s |

No regression. The second enzyme costs ~50%, split between the extra search and
a 2.5× larger output to write.

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
