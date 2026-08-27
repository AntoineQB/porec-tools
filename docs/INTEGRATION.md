# Using `porec` inside wf-pore-c

`wf-pore-c` shells out to `pore-c-py digest "$cutter"` in its
`digest_align_annotate` process. Making it multi-enzyme is a one-line change
plus getting the package into the container.

Three routes, from least to most invasive.

---

## A. Run the digest yourself, keep the workflow for the rest

The simplest option, and the one to prefer for a one-off reanalysis: digest
with `porec`, then feed the monomers to the rest of your pipeline.

```bash
porec digest DpnII,NlaIII reads.bam \
    --output monomers.ns.bam \
    --stats  digest_stats.tsv \
    --threads 8
```

Add `--max_monomers 250 --excluded_list dropped.txt` to mirror what the
workflow does with over-digested reads.

Nothing in the workflow is touched, and `digest_stats.tsv` records how many
sites each enzyme had and how much of the fragmentation it accounts for. (It
does *not* tell you whether an enzyme cut, see the last section.)

---

## B. Patch the workflow

`patches/wf-pore-c-multicutter.patch` applies to the workflow checkout:

```bash
nextflow pull epi2me-labs/wf-pore-c
cd ~/.nextflow/assets/epi2me-labs/wf-pore-c
git apply /path/to/porec-tools/patches/wf-pore-c-multicutter.patch
```

It does two things:

1. swaps `pore-c-py digest` for `porec digest` in `modules/local/pore-c.nf`
   at both call sites;
2. documents the comma-separated form in `nextflow_schema.json`.

The two call sites matter. `digest_align_annotate` branches on `chunk_size`,
and the two branches order the positionals differently:

```bash
# chunked: enzyme only, BAM arrives on stdin
bamindex fetch --chunk=N concatemers.bam | pore-c-py digest "${meta.cutter}" ...
# not chunked: input file first, enzyme second
pore-c-py digest "concatemers.bam" "${meta.cutter}" ...
```

`porec` accepts either order, and implements the `--max_monomers` and
`--excluded_list` options both branches pass. Patching only one branch, or
dropping those options, leaves the workflow broken on one of its two paths.

The container must then contain the package, see below.

---

## C. Build a container with the package installed

```dockerfile
FROM ontresearch/wf-pore-c:sha3787c234c0cacf66a67fb77da223cc2e1cb0baf0
USER root
RUN pip install --no-cache-dir git+https://github.com/YOUR-USERNAME/porec-tools.git
USER $WF_UID
```

```bash
docker build -t wf-porec:1.0.0 .
nextflow run epi2me-labs/wf-pore-c \
    --cutter 'DpnII,NlaIII' \
    -process.container wf-porec:1.0.0 \
    ...
```

`pore-c-py` stays installed and untouched, so you can compare the two on the
same input.

---

## Checking it worked

Before committing to a long run, resolve the enzymes without reading data:

```bash
porec digest DpnII,NlaIII --dry-run
```

```
Enzymes resolved from 'DpnII,NlaIII':
  enzyme  site  cut (both strands)  sticky end  1 site per  recognition
  DpnII   GATC  N^GATC_N            5' GATC     256 bp      palindromic
  NlaIII  CATG  _CATG^              3' CATG     256 bp      palindromic
```

Check the site and the cut position against your protocol here. It costs a
second, and a wrong enzyme wastes the whole run.

After the digest, read the per-enzyme table. `0 sites` for an enzyme means
the name is wrong or the reads are not from that run, and is flagged as a
warning. A large count means only that the motif is present, which it always
is. See the caveat below.

Provenance is recorded in the BAM header:

```bash
samtools view -H monomers.ns.bam | grep porec
# @PG ID:porec PN:porec VN:1.0.0 CL:... DS:multi-enzyme digest: DpnII(GATC), NlaIII(CATG)
```

---

## After the workflow: undo the cuts the enzyme never made

The digest cuts at every recognition site on purpose,
because from the read sequence alone an uncut site inside a fragment is
indistinguishable from a reconstituted ligation junction. That over-cutting has
to be undone once the monomers are aligned, or the Hi-C diagonal is inflated by
contacts that never happened.

```bash
porec merge sample.ns.bam \
    --output fragments.tsv.gz \
    --pairs contacts.pairs --sizes hg38.sizes.genome \
    --stats merge_stats.json --min-fragments 2
```

Where it sits in the chain:

```
concatemers.bam
  -> porec digest      cut at every site (over-cuts, unavoidably)
  -> minimap2               each monomer gets a genomic position
  -> porec merge       glue back what was never really cut
  -> .pairs -> cooler / juicer
```

`fragments.tsv.gz` is one row per merged fragment: read, chromosome, start,
end, midpoint, strand, MAPQ, and how many monomers were glued. `contacts.pairs`
is 4DN v1.0, ready for `cooler cload pairs` or `juicer_tools pre`.

Read the before/after table it prints. A large `<1kb` share before merging
is the artefact itself: pieces of one uncut restriction fragment being paired
with each other, landing on the diagonal.

---

## Which enzymes did my library actually use?

The digest report cannot answer this, and does not pretend to: it counts sites
present in the reads, and every motif is present by chance. The test that does
work needs aligned monomers, and ships with this package:

```bash
porec-junctions monomers.aligned.ns.bam hg38.fa \
    --enzymes DpnII,NlaIII,HindIII,HinfI
```

It takes the boundaries between consecutive monomers that jump more than 1 kb,
which are genuine ligation junctions, and asks how often each enzyme's site sits there,
against random positions on the same chromosomes:

```
30,000 junction boundaries (MAPQ >= 20, jump > 1,000 bp)

  enzyme     site      at junctions   at random  enrichment   verdict
  DpnII      GATC             39.5%        1.3%       29.5x   yes, main enzyme
  NlaIII     CATG             10.3%        2.3%        4.4x   yes, secondary
  HindIII    AAGCTT            0.0%        0.1%        0.1x   no
  HinfI      GANTC             1.3%        1.5%        0.9x   no
```

That is the run this tool was written for. The protocol notes said HindIII;
`AAGCTT` sits below background at junctions, so HindIII never cut. The
second enzyme was NlaIII. Same verdicts on the second sample of the pair.

Notes on reading it:

- Use the enrichment column, not the percentages. Aligners soft-clip a few
  bases at monomer ends, so the boundary drifts off the true cut and the
  absolute rates understate. `--tol 10` recovers DpnII to ~79% but inflates the
  random background too, dropping every enrichment.
- The outer ends of a concatemer are ignored. They carry adapters, not cut
  sites. (A diagnostic based on read termini was tried first and returns 0% for
  every enzyme, including ones that certainly cut, for exactly this reason.)
- Running this before the digest saves re-running it. It only needs an
  aligned BAM from any earlier single-enzyme run.
