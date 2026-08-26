# Using `pore-c-aqb` inside wf-pore-c

`wf-pore-c` shells out to `pore-c-py digest "$cutter"` in its
`digest_align_annotate` process. Making it multi-enzyme is a one-line change
plus getting the package into the container.

Three routes, from least to most invasive.

---

## A. Run the digest yourself, keep the workflow for the rest

The simplest option, and the one to prefer for a one-off reanalysis: digest
with `pore-c-aqb`, then feed the monomers to the rest of your pipeline.

```bash
pore-c-aqb digest DpnII,NlaIII reads.bam \
    --output monomers.ns.bam \
    --stats  digest_stats.tsv \
    --threads 8
```

Nothing in the workflow is touched, and `digest_stats.tsv` tells you straight
away whether both enzymes actually cut.

---

## B. Patch the workflow

`patches/wf-pore-c-multicutter.patch` applies to the workflow checkout:

```bash
nextflow pull epi2me-labs/wf-pore-c
cd ~/.nextflow/assets/epi2me-labs/wf-pore-c
git apply /path/to/wf-pore-c_AQB/patches/wf-pore-c-multicutter.patch
```

It does two things:

1. swaps `pore-c-py digest` for `pore-c-aqb digest` in `modules/local/pore-c.nf`;
2. drops the single-value validation on `--cutter` in `nextflow_schema.json`,
   so `--cutter 'DpnII,NlaIII'` is accepted.

The container must then contain the package — see below.

---

## C. Build a container with the package installed

```dockerfile
FROM ontresearch/wf-pore-c:sha3787c234c0cacf66a67fb77da223cc2e1cb0baf0
USER root
RUN pip install --no-cache-dir git+https://github.com/AQB/wf-pore-c_AQB.git
USER $WF_UID
```

```bash
docker build -t wf-pore-c-aqb:1.0.0 .
nextflow run epi2me-labs/wf-pore-c \
    --cutter 'DpnII,NlaIII' \
    -process.container wf-pore-c-aqb:1.0.0 \
    ...
```

`pore-c-py` stays installed and untouched, so you can compare the two on the
same input.

---

## Checking it worked

**Before committing to a long run**, resolve the enzymes without reading data:

```bash
pore-c-aqb digest DpnII,NlaIII --dry-run
# DpnII    GATC    fst5=0    palindromic=True
# NlaIII   CATG    fst5=4    palindromic=True
```

**After the digest**, read the per-enzyme report. An enzyme at `0 cuts` either
was not in the reaction, failed, or is a typo:

```
  DpnII        GATC             61,229 cuts ( 39.1%)
  NlaIII       CATG             95,424 cuts ( 61.0%)
```

**Provenance** is recorded in the BAM header:

```bash
samtools view -H monomers.ns.bam | grep pore-c-aqb
# @PG ID:pore-c-aqb PN:pore-c-aqb VN:1.0.0 CL:... DS:multi-enzyme digest: DpnII(GATC), NlaIII(CATG)
```

---

## Which enzymes did my library actually use?

If the protocol is unclear, the data can be asked directly. Take the junctions
between consecutive aligned monomers that jump more than 1 kb, and look at what
motif sits at the boundary:

```bash
samtools view monomers.ns.bam | python3 - <<'PY'
# see the project this tool came from for a fuller version;
# the principle: extract junction boundaries, then compare motif frequency
# at those positions against random positions on the same chromosomes.
PY
```

On the dataset this tool was written for, that analysis gave `GATC` at 86% of
junctions (27× above random) and `CATG` at 21% (16× above random), while
`AAGCTT` sat at 0.1% — below background. The protocol notes said HindIII; the
data said NlaIII.

Doing this **before** the digest saves re-running it.
