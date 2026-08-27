# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""A progress bar, so a long run tells you where it is and when it will end.

Digesting a full library takes tens of minutes and merging a 24 GB BAM takes a
couple of them. With no output during that time you cannot tell a slow run
from a stuck one, and the only options are to wait or to kill it.

How the percentage is real, not guessed
---------------------------------------
The number of records in a BAM is not known without reading it, so counting
records cannot give a percentage. But a BAM is BGZF-compressed, and
``AlignmentFile.tell()`` returns a virtual offset whose top 48 bits are the
position in the compressed file. Compared with the file size, that is a genuine
fraction of the work done, and the ETA follows from it.

Reading from a pipe has no size, so the bar falls back to counting and rate
only - honest about not knowing rather than inventing a total.

Deliberate choices
------------------
* stderr, never stdout. Some commands write their real output to stdout;
  a bar there would corrupt it.
* Off unless stderr is a terminal. Redirect to a file or a Nextflow log and
  you get no carriage-return soup. ``--progress`` forces it on anyway.
* No dependency. tqdm would do this well, but adding a dependency to a tool
  people install into an existing pipeline is a cost, for forty lines.
* Rate-limited redraws. At most every 0.15 s, so the bar never becomes the
  bottleneck it is reporting on.
"""
from __future__ import annotations

import os
import shutil
import sys
import time

__all__ = ["Progress", "add_progress_arguments", "progress_enabled"]

_MIN_REDRAW_INTERVAL = 0.15
_BAR_WIDTH = 24


def add_progress_arguments(parser) -> None:
    """``--progress`` / ``--no-progress``, with the same meaning everywhere."""
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--progress", dest="progress", action="store_const", const=True,
        default=None,
        help=("Show a progress bar on stderr. On by default when stderr is a "
              "terminal; use this to force it when redirecting to a file."))
    group.add_argument(
        "--no-progress", dest="progress", action="store_const", const=False,
        help="Never show a progress bar.")


def progress_enabled(args) -> bool:
    """Whether to draw: --progress/--no-progress first, else "is this a TTY?".

    ``--quiet`` deliberately does *not* turn the bar off. It controls logging,
    and a bar is not a log line: it exists only on a terminal, where somebody
    is watching and wants to know how long this will take. Silencing both from
    one flag also made the two commands disagree, since only ``digest`` has a
    log level. ``--no-progress`` is the way to turn it off.
    """
    choice = getattr(args, "progress", None)
    if choice is not None:
        return choice
    return sys.stderr.isatty()


def _format_duration(seconds: float) -> str:
    if seconds != seconds or seconds in (float("inf"), float("-inf")):
        return "--:--"
    seconds = int(max(0, seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes}:{secs:02d}"


class Progress:
    """Draws a single self-updating line while work happens.

    Usable as a context manager; the line is cleared on the way out so it never
    interleaves with the report that follows.

    :param label: what is happening, e.g. ``"digesting"``
    :param unit: what is being counted, e.g. ``"reads"``
    :param total: total units, when known
    :param position: callable returning bytes consumed so far
    :param size: total bytes, when ``position`` is given
    :param enabled: draw at all
    """

    def __init__(self, label: str, unit: str = "records", total=None,
                 position=None, size=None, enabled: bool = True,
                 stream=None):
        self.label = label
        self.unit = unit
        self.total = total
        self.position = position
        self.size = size if size and size > 0 else None
        self.stream = stream if stream is not None else sys.stderr
        self.enabled = bool(enabled)
        self.count = 0
        self._started = time.monotonic()
        self._last_draw = 0.0
        self._drawn = False

    # -- context manager ---------------------------------------------------
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    # -- driving it --------------------------------------------------------
    def update(self, n: int = 1) -> None:
        self.count += n
        if not self.enabled:
            return
        now = time.monotonic()
        if now - self._last_draw < _MIN_REDRAW_INTERVAL:
            return
        self._last_draw = now
        self._draw(now)

    def close(self, summary: str | None = None) -> None:
        """Clear the line, optionally leaving one final summary behind."""
        if not self.enabled:
            return
        if self._drawn:
            self.stream.write("\r" + " " * self._width() + "\r")
        if summary:
            self.stream.write(summary.rstrip() + "\n")
        self.stream.flush()
        self._drawn = False

    # -- drawing -----------------------------------------------------------
    def _width(self) -> int:
        return max(40, min(shutil.get_terminal_size((100, 24)).columns, 120))

    def fraction(self):
        """How far along, in [0, 1], or None when that cannot be known."""
        if self.total:
            return min(1.0, self.count / self.total)
        if self.position is not None and self.size:
            try:
                return min(1.0, max(0.0, self.position() / self.size))
            except Exception:            # pragma: no cover - defensive
                return None
        return None

    def render(self, now: float | None = None) -> str:
        now = time.monotonic() if now is None else now
        elapsed = now - self._started
        rate = self.count / elapsed if elapsed > 0 else 0.0
        done = self.fraction()

        parts = [self.label]
        if done is not None:
            filled = int(round(done * _BAR_WIDTH))
            parts.append("[" + "#" * filled + "." * (_BAR_WIDTH - filled) + "]")
            parts.append(f"{done * 100:3.0f}%")
        parts.append(f"{self.count:,} {self.unit}")
        parts.append(f"{rate:,.0f}/s")
        parts.append(f"elapsed {_format_duration(elapsed)}")
        if done is not None and 0 < done < 1 and rate > 0:
            parts.append(f"eta {_format_duration(elapsed * (1 - done) / done)}")
        return "  ".join(parts)

    def _draw(self, now: float) -> None:
        line = self.render(now)
        width = self._width()
        if len(line) > width:
            line = line[:width]
        self.stream.write("\r" + line.ljust(width))
        self.stream.flush()
        self._drawn = True


def bam_progress(label: str, path: str, handle, args, unit: str = "records"):
    """A :class:`Progress` tracking how far through a BAM file we are.

    Falls back to counting without a percentage when the input is a pipe, which
    has no size to measure against.
    """
    size = None
    if path and path != "-" and os.path.exists(path):
        size = os.path.getsize(path)
    return Progress(
        label, unit=unit,
        position=(lambda: handle.tell() >> 16) if size else None,
        size=size,
        enabled=progress_enabled(args),
    )
