"""Command-line behaviour, including drop-in compatibility."""
from __future__ import annotations

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


def test_stats_file_reports_each_enzyme(tmp_path, concatemers):
    out, stats = tmp_path / "out.bam", tmp_path / "stats.tsv"
    main(["digest", "DpnII,NlaIII", concatemers, "--output", str(out),
          "--stats", str(stats), "--quiet"])
    text = stats.read_text()
    assert "cuts_DpnII" in text and "cuts_NlaIII" in text
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
