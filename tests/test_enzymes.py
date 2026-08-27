"""Enzyme parsing, resolution and validation."""
import pytest

from porec_tools.enzymes import (
    EnzymeSpecError,
    describe_enzymes,
    parse_enzyme_spec,
    resolve_enzymes,
)


@pytest.mark.parametrize("spec,expected", [
    ("DpnII", ["DpnII"]),
    ("DpnII,NlaIII", ["DpnII", "NlaIII"]),
    ("DpnII+NlaIII", ["DpnII", "NlaIII"]),
    ("DpnII NlaIII", ["DpnII", "NlaIII"]),
    ("DpnII;NlaIII", ["DpnII", "NlaIII"]),
    ("  DpnII , NlaIII ", ["DpnII", "NlaIII"]),
    (["DpnII", "NlaIII"], ["DpnII", "NlaIII"]),
])
def test_parse_accepts_common_separators(spec, expected):
    assert parse_enzyme_spec(spec) == expected


def test_parse_preserves_order_and_drops_duplicates():
    assert parse_enzyme_spec("NlaIII,DpnII,NlaIII") == ["NlaIII", "DpnII"]


@pytest.mark.parametrize("spec", ["", "   ", ",,,", []])
def test_parse_rejects_empty(spec):
    with pytest.raises(EnzymeSpecError):
        parse_enzyme_spec(spec)


def test_resolve_single():
    (enz,) = resolve_enzymes("DpnII")
    assert enz.name == "DpnII"
    assert enz.site == "GATC"
    assert enz.fst5 == 0


def test_resolve_multiple_keeps_order():
    got = resolve_enzymes("DpnII,NlaIII")
    assert [e.name for e in got] == ["DpnII", "NlaIII"]
    assert [e.site for e in got] == ["GATC", "CATG"]


def test_resolve_collapses_isoschizomers():
    """DpnII and MboI recognise GATC and cut identically: one is enough."""
    got = resolve_enzymes("DpnII,MboI")
    assert len(got) == 1, "isoschizomers should not be searched twice"
    assert got[0].site == "GATC"


def test_unknown_enzyme_is_rejected_with_a_suggestion():
    with pytest.raises(EnzymeSpecError) as exc:
        resolve_enzymes("DpnI I")
    assert "Unknown restriction enzyme" in str(exc.value)


def test_case_mistake_suggests_the_right_name():
    with pytest.raises(EnzymeSpecError) as exc:
        resolve_enzymes("dpnii")
    assert "DpnII" in str(exc.value)


def test_all_names_validated_before_use():
    """A typo in the second enzyme must fail immediately, not mid-run."""
    with pytest.raises(EnzymeSpecError) as exc:
        resolve_enzymes("DpnII,NotAnEnzyme")
    assert "NotAnEnzyme" in str(exc.value)


def test_describe():
    assert describe_enzymes(resolve_enzymes("DpnII,NlaIII")) == \
        "DpnII(GATC), NlaIII(CATG)"


def test_ambiguity_codes_supported():
    """HinfI is GANTC: an IUPAC-containing site must resolve."""
    (enz,) = resolve_enzymes("HinfI")
    assert enz.site == "GANTC"
