"""Command-line behaviour, including drop-in compatibility."""
from __future__ import annotations

import subprocess
import sys

import pysam
import pytest

from pore_c_aqb.cli import main

from conftest import make_read, random_seq, write_bam


@pytest.fixture
def concatemers(tmp_path):
    reads = [make_read(f"read{i}", random_seq(500 + 97 * i, seed=i))
             for i in range(10)]
    return write_bam(tmp_path / "in.bam", reads)


def _monomers(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as fh:
        return list(fh)


def test_single_enzyme_runs(tmp_path, concatemers):
    out = tmp_path / "out.bam"
    assert main(["digest", "DpnII", concatemers, "--output", str(out),
                 "--quiet"]) == 0
    assert len(_monomers(out)) > 10


def test_two_enzymes_give_at_least_as_many_monomers(tmp_path, concatemers):
    one, two = tmp_path / "one.bam", tmp_path / "two.bam"
    main(["digest", "DpnII", concatemers, "--output", str(one), "--quiet"])
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(two),
          "--quiet"])
    assert len(_monomers(two)) >= len(_monomers(one))


def test_enzyme_order_does_not_change_the_result(tmp_path, concatemers):
    """The digest is a set union: order is irrelevant."""
    a, b = tmp_path / "a.bam", tmp_path / "b.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(a), "--quiet"])
    main(["digest", "NlaIII,DpnII", concatemers, "--output", str(b), "--quiet"])
    names_a = [r.query_name for r in _monomers(a)]
    names_b = [r.query_name for r in _monomers(b)]
    assert names_a == names_b


def test_unknown_enzyme_exits_with_code_2(tmp_path, concatemers, caplog):
    out = tmp_path / "out.bam"
    assert main(["digest", "NotAnEnzyme", concatemers, "--output", str(out),
                 "--quiet"]) == 2
    assert not out.exists(), "no output should be written on a bad enzyme"


def test_dry_run_lists_enzymes(capsys):
    assert main(["digest", "DpnII,NlaIII", "--dry-run", "--quiet"]) == 0
    out = capsys.readouterr().out
    assert "DpnII" in out and "GATC" in out
    assert "NlaIII" in out and "CATG" in out
    # the cut must be shown as a readable position, not as "fst5=0"
    assert "N^GATC_N" in out and "_CATG^" in out
    assert "fst5" not in out


def test_stats_file_reports_each_enzyme(tmp_path, concatemers):
    out, stats = tmp_path / "out.bam", tmp_path / "stats.tsv"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(out),
          "--stats", str(stats), "--quiet"])
    text = stats.read_text()
    assert "sites_found_DpnII" in text and "sites_found_NlaIII" in text
    assert "cuts_DpnII" not in text, "the TSV must not call site counts cuts"
    assert "chance_spacing_bp_DpnII" in text
    assert "observed_spacing_bp_DpnII" in text
    assert "overhang_DpnII" in text
    assert "concatemers\t10" in text


def test_provenance_recorded_in_bam_header(tmp_path, concatemers):
    out = tmp_path / "out.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(out),
          "--quiet"])
    with pysam.AlignmentFile(str(out), "rb", check_sq=False) as fh:
        pg = fh.header.to_dict().get("PG", [])
    entry = [p for p in pg if p.get("ID") == "pore-c-aqb"]
    assert entry, "a @PG line must record how the file was made"
    assert "DpnII" in entry[0]["DS"] and "NlaIII" in entry[0]["DS"]


def test_monomers_reconstruct_the_input(tmp_path, concatemers):
    """Nothing is lost or duplicated: monomers concatenate back to the reads."""
    out = tmp_path / "out.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(out),
          "--quiet"])
    rebuilt: dict[str, str] = {}
    for r in _monomers(out):
        rebuilt.setdefault(r.get_tag("MI"), "")
        rebuilt[r.get_tag("MI")] += r.query_sequence
    with pysam.AlignmentFile(concatemers, "rb", check_sq=False) as fh:
        for read in fh:
            assert rebuilt[read.query_name] == read.query_sequence


def test_stdin_without_header_is_refused(tmp_path):
    with pytest.raises(SystemExit):
        main(["digest", "DpnII", "-", "--output", str(tmp_path / "o.bam"),
              "--quiet"])


# --------------------------------------------------------------------------
# compatibility with how wf-pore-c actually invokes the digest
# --------------------------------------------------------------------------
def test_accepts_input_before_enzyme(tmp_path, concatemers):
    """wf-pore-c's non-chunked branch calls `digest concatemers.bam "$cutter"`.

    Its chunked branch calls `digest "$cutter"` on stdin. Both orders have to
    work or the patched workflow fails on one of its two code paths.
    """
    a, b = tmp_path / "a.bam", tmp_path / "b.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(a), "--quiet"])
    main(["digest", concatemers, "DpnII,NlaIII", "--output", str(b), "--quiet"])
    with pysam.AlignmentFile(str(a), check_sq=False) as fa, \
         pysam.AlignmentFile(str(b), check_sq=False) as fb:
        assert [r.query_name for r in fa] == [r.query_name for r in fb]


def test_max_reads_takes_exactly_n(tmp_path, concatemers, capsys):
    out = tmp_path / "out.bam"
    main(["digest", "DpnII", concatemers, "--output", str(out),
          "--max_reads", "3", "--quiet"])
    with pysam.AlignmentFile(str(out), check_sq=False) as f:
        assert len({r.get_tag("MI") for r in f}) == 3


def test_max_monomers_excludes_whole_reads(tmp_path, concatemers):
    """An over-digested read is dropped entire, not trimmed."""
    out, lst = tmp_path / "out.bam", tmp_path / "excluded.txt"
    excl = tmp_path / "excluded.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(out),
          "--max_monomers", "2", "--excluded_list", str(lst),
          "--excluded_bam", str(excl), "--quiet"])
    names = lst.read_text().split()
    assert names, "the fixture should contain reads with more than 2 monomers"
    with pysam.AlignmentFile(str(excl), check_sq=False) as f:
        assert sorted(r.query_name for r in f) == sorted(names)
    with pysam.AlignmentFile(str(out), check_sq=False) as f:
        kept = {r.get_tag("MI") for r in f}
    assert kept.isdisjoint(names), "an excluded read must emit no monomer"


def test_excluded_reads_are_whole_and_untrimmed(tmp_path, concatemers):
    excl = tmp_path / "excluded.bam"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(tmp_path / "o.bam"),
          "--max_monomers", "1", "--excluded_bam", str(excl), "--quiet"])
    with pysam.AlignmentFile(concatemers, check_sq=False) as f:
        original = {r.query_name: r.query_sequence for r in f}
    with pysam.AlignmentFile(str(excl), check_sq=False) as f:
        for r in f:
            assert r.query_sequence == original[r.query_name]


def test_directory_input_is_expanded(tmp_path, concatemers):
    import shutil
    d = tmp_path / "in"
    d.mkdir()
    shutil.copy(concatemers, d / "one.bam")
    shutil.copy(concatemers, d / "two.bam")
    out = tmp_path / "out.bam"
    main(["digest", "DpnII", str(d), "--glob", "*.bam",
          "--output", str(out), "--quiet"])
    with pysam.AlignmentFile(str(out), check_sq=False) as f:
        n_dir = sum(1 for _ in f)
    single = tmp_path / "single.bam"
    main(["digest", "DpnII", concatemers, "--output", str(single), "--quiet"])
    with pysam.AlignmentFile(str(single), check_sq=False) as f:
        n_one = sum(1 for _ in f)
    assert n_dir == 2 * n_one


def test_empty_directory_is_reported(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    with pytest.raises(SystemExit, match="No file matching"):
        main(["digest", "DpnII", str(d), "--output", str(tmp_path / "o.bam")])


def test_missing_input_is_reported(tmp_path):
    with pytest.raises(SystemExit, match="Input file not found"):
        main(["digest", "DpnII", str(tmp_path / "nope.bam"),
              "--output", str(tmp_path / "o.bam")])


def test_enzyme_named_like_a_file_is_still_an_enzyme(tmp_path, concatemers):
    """The swap only triggers when the first positional really is a file."""
    out = tmp_path / "o.bam"
    assert main(["digest", "DpnII", concatemers,
                 "--output", str(out), "--quiet"]) == 0


def _cli(*argv, stdin=None):
    """Run the CLI in a subprocess, so stdin is a real pipe."""
    return subprocess.run(
        [sys.executable, "-m", "pore_c_aqb.cli", *argv],
        input=stdin, capture_output=True)


def test_digest_reads_from_a_real_stdin_pipe(tmp_path, concatemers):
    """wf-pore-c's chunked branch pipes the BAM in on stdin.

    A pipe can only be opened once, which an in-process test using a file path
    never exercises: the header has to come from the same handle the reads do.
    """
    out = tmp_path / "out.bam"
    with open(concatemers, "rb") as fh:
        data = fh.read()
    proc = _cli("digest", "DpnII,NlaIII", "-", "--header", concatemers,
                "--output", str(out), "--quiet", stdin=data)
    assert proc.returncode == 0, proc.stderr.decode()[-800:]
    with pysam.AlignmentFile(str(out), check_sq=False) as f:
        from_pipe = [r.query_name for r in f]

    direct = tmp_path / "direct.bam"
    main(["digest", "DpnII,NlaIII", concatemers,
          "--output", str(direct), "--quiet"])
    with pysam.AlignmentFile(str(direct), check_sq=False) as f:
        assert from_pipe == [r.query_name for r in f]


def test_stdin_works_with_the_workflow_argument_order(tmp_path, concatemers):
    """`digest "$cutter"` with no input positional, exactly as wf-pore-c does."""
    out = tmp_path / "out.bam"
    with open(concatemers, "rb") as fh:
        data = fh.read()
    proc = _cli("digest", "DpnII,NlaIII", "--header", concatemers,
                "--max_monomers", "250",
                "--excluded_list", str(tmp_path / "filtered_reads.txt"),
                "--output", str(out), "--quiet", stdin=data)
    assert proc.returncode == 0, proc.stderr.decode()[-800:]
    with pysam.AlignmentFile(str(out), check_sq=False) as f:
        assert sum(1 for _ in f) > 0
