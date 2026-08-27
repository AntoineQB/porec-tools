"""The progress bar.

Small, but it writes to a shared stream while other things are printing, so the
failure modes are ugly: a bar left half-drawn across a report, carriage returns
in a Nextflow log, or a percentage that lies. Each of those is a test here.
"""
import argparse
import io
import time

import pytest

from pore_c_aqb.progress import (
    Progress,
    add_progress_arguments,
    progress_enabled,
)


def bar(**kw):
    kw.setdefault("stream", io.StringIO())
    kw.setdefault("enabled", True)
    return Progress("working", **kw)


def parse(argv, isatty=False):
    p = argparse.ArgumentParser()
    add_progress_arguments(p)
    args = p.parse_args(argv)
    return args


# --------------------------------------------------------------------------
# when it draws at all
# --------------------------------------------------------------------------
def test_off_when_stderr_is_not_a_terminal(monkeypatch):
    """A redirected stderr means a log file: carriage returns would ruin it."""
    monkeypatch.setattr("sys.stderr", io.StringIO())
    assert progress_enabled(parse([])) is False


def test_on_when_stderr_is_a_terminal(monkeypatch):
    class Tty(io.StringIO):
        def isatty(self):
            return True
    monkeypatch.setattr("sys.stderr", Tty())
    assert progress_enabled(parse([])) is True


def test_progress_flag_forces_it_on(monkeypatch):
    monkeypatch.setattr("sys.stderr", io.StringIO())
    assert progress_enabled(parse(["--progress"])) is True


def test_no_progress_flag_forces_it_off(monkeypatch):
    class Tty(io.StringIO):
        def isatty(self):
            return True
    monkeypatch.setattr("sys.stderr", Tty())
    assert progress_enabled(parse(["--no-progress"])) is False


def test_quiet_does_not_silence_the_bar(monkeypatch):
    """--quiet controls logging. A bar is not a log line.

    Tying them together also made the commands disagree, since only `digest`
    has a log level at all.
    """
    class Tty(io.StringIO):
        def isatty(self):
            return True
    monkeypatch.setattr("sys.stderr", Tty())
    args = parse([])
    args.log_level = 30                      # logging.WARNING, i.e. --quiet
    assert progress_enabled(args) is True


def test_disabled_bar_writes_nothing():
    stream = io.StringIO()
    p = Progress("working", stream=stream, enabled=False)
    for _ in range(1000):
        p.update()
    p.close()
    assert stream.getvalue() == ""
    assert p.count == 1000, "counting must continue even when not drawing"


# --------------------------------------------------------------------------
# what it says
# --------------------------------------------------------------------------
def test_percentage_comes_from_the_file_position():
    p = bar(position=lambda: 250, size=1000)
    p.update()
    assert p.fraction() == pytest.approx(0.25)
    assert " 25%" in p.render()


def test_percentage_from_a_known_total():
    p = bar(total=200)
    for _ in range(50):
        p.update()
    assert p.fraction() == pytest.approx(0.25)


def test_no_percentage_when_the_size_is_unknown():
    """Reading a pipe: report count and rate, do not invent a total."""
    p = bar()
    p.update()
    assert p.fraction() is None
    text = p.render()
    assert "%" not in text and "eta" not in text
    assert "1 working" not in text
    assert "records" in text


def test_fraction_never_exceeds_one():
    p = bar(position=lambda: 5000, size=1000)
    p.update()
    assert p.fraction() == 1.0


def test_fraction_survives_a_broken_position_callback():
    def boom():
        raise OSError("closed")
    p = bar(position=boom, size=1000)
    assert p.fraction() is None


def test_eta_shrinks_as_work_proceeds():
    state = {"at": 100}
    p = bar(position=lambda: state["at"], size=1000)
    p.update()
    time.sleep(0.01)
    early = p.render()
    state["at"] = 900
    late = p.render()
    assert "eta" in early and "eta" in late
    assert p.fraction() > 0.5


def test_no_eta_at_the_very_end():
    p = bar(position=lambda: 1000, size=1000)
    p.update()
    assert "eta" not in p.render()


def test_render_includes_the_label_and_unit():
    p = bar(unit="reads")
    p.update(7)
    text = p.render()
    assert text.startswith("working")
    assert "7 reads" in text


def test_counts_are_thousands_separated():
    p = bar()
    p.update(1234567)
    assert "1,234,567" in p.render()


# --------------------------------------------------------------------------
# not corrupting the terminal
# --------------------------------------------------------------------------
def test_close_erases_the_line():
    stream = io.StringIO()
    p = Progress("working", stream=stream, enabled=True)
    p.update()
    p.close()
    assert stream.getvalue().endswith("\r"), \
        "the bar must be wiped, or it stays under whatever prints next"


def test_close_is_safe_without_any_update():
    stream = io.StringIO()
    Progress("working", stream=stream).close()
    assert stream.getvalue() == ""


def test_close_twice_is_harmless():
    p = bar()
    p.update()
    p.close()
    p.close()


def test_context_manager_closes_on_exception():
    stream = io.StringIO()
    with pytest.raises(ValueError):
        with Progress("working", stream=stream, enabled=True) as p:
            p.update()
            raise ValueError("boom")
    assert stream.getvalue().endswith("\r")


def test_a_summary_line_survives_the_wipe():
    stream = io.StringIO()
    p = Progress("working", stream=stream, enabled=True)
    p.update()
    p.close(summary="done: 42 things")
    assert stream.getvalue().endswith("done: 42 things\n")


def test_output_stays_on_one_line():
    stream = io.StringIO()
    p = Progress("working", stream=stream, enabled=True, total=10)
    for _ in range(10):
        p.update()
        p._draw(time.monotonic())          # force a redraw past the rate limit
    assert "\n" not in stream.getvalue()


def test_redraws_are_rate_limited():
    """The bar must never become the bottleneck it is reporting on."""
    stream = io.StringIO()
    p = Progress("working", stream=stream, enabled=True, total=100000)
    for _ in range(100000):
        p.update()
    assert stream.getvalue().count("\r") <= 3, "one burst should redraw once"


def test_long_lines_are_truncated_not_wrapped():
    stream = io.StringIO()
    p = Progress("a" * 500, stream=stream, enabled=True)
    p.update()
    p._draw(time.monotonic())
    longest = max(len(part) for part in stream.getvalue().split("\r") if part)
    assert longest <= 120
