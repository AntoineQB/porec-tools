# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Command-line interface, argument-compatible with ``pore-c-py digest``.

The only difference is that the positional ``enzyme`` argument accepts several
names, so an existing wf-pore-c invocation keeps working unchanged:

    pore-c-aqb digest DpnII        input.bam --output monomers.bam
    pore-c-aqb digest DpnII,NlaIII input.bam --output monomers.bam
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pysam

from pore_c_aqb import __version__
from pore_c_aqb.digest import DigestStats, get_concatemer_seqs
from pore_c_aqb.enzymes import (
    EnzymeSpecError,
    describe_enzymes,
    resolve_enzymes,
)

PROG = "pore-c-aqb"
logger = logging.getLogger(PROG)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description=(
            "Digest concatemers into monomers using one or more restriction "
            "enzymes. Drop-in replacement for 'pore-c-py digest' that accepts "
            "a comma-separated list of enzymes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"{PROG} {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser(
        "digest",
        help="Digest concatemers into monomers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    d.add_argument(
        "enzyme",
        help=("One or more restriction enzyme names. Separate several with "
              "commas, e.g. 'DpnII,NlaIII'. Names are those of Biopython's "
              "Restriction module."),
    )
    d.add_argument(
        "input", nargs="?", default="-",
        help="Unaligned BAM of concatemers ('-' for stdin).")
    d.add_argument(
        "--output", default="-",
        help="Output unaligned BAM ('-' for stdout).")
    d.add_argument(
        "--header", type=Path, default=None,
        help=("BAM whose header should be copied to the output. Required when "
              "reading from stdin, which carries no accessible header."))
    d.add_argument(
        "--remove_tags", nargs="+", default=None,
        help="Additional SAM tags to strip from the output.")
    d.add_argument(
        "--max_reads", type=int, default=0,
        help="Stop after this many concatemers (0 = all). Useful for testing.")
    d.add_argument(
        "--threads", type=int, default=1,
        help="Compute threads for BAM compression.")
    d.add_argument(
        "--stats", type=Path, default=None,
        help="Write the per-enzyme cut report to this file as TSV.")
    d.add_argument(
        "--dry-run", action="store_true",
        help="Resolve and print the enzymes, then exit without reading input.")

    verb = d.add_mutually_exclusive_group()
    verb.add_argument(
        "--debug", action="store_const", dest="log_level",
        const=logging.DEBUG, default=logging.INFO,
        help="Verbose logging of debug information.")
    verb.add_argument(
        "--quiet", action="store_const", dest="log_level",
        const=logging.WARNING, default=logging.INFO,
        help="Minimal logging; warnings only.")
    d.add_argument("--logfile", type=Path, default=None,
                   help="Write logs to this file as well as stderr.")
    return parser


def _setup_logging(level: int, logfile: Path | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile:
        handlers.append(logging.FileHandler(logfile))
    logging.basicConfig(
        level=level, handlers=handlers,
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def _open_input(path: str, header_path: Path | None):
    if path == "-":
        if header_path is None:
            raise SystemExit(
                "Reading from stdin requires --header: a piped BAM carries no "
                "header that can be copied to the output."
            )
        return pysam.AlignmentFile("-", "rb", check_sq=False)
    return pysam.AlignmentFile(path, "rb", check_sq=False)


def _output_header(infile, header_path: Path | None, enzymes) -> dict:
    """Header for the output, with a provenance line for the digest."""
    if header_path is not None:
        with pysam.AlignmentFile(str(header_path), "rb", check_sq=False) as h:
            header = h.header.to_dict()
    else:
        header = infile.header.to_dict()

    header.setdefault("PG", []).append({
        "ID": PROG,
        "PN": PROG,
        "VN": __version__,
        "CL": " ".join(sys.argv),
        "DS": f"multi-enzyme digest: {describe_enzymes(enzymes)}",
    })
    return header


def _write_stats(path: Path, enzymes, stats: DigestStats) -> None:
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"enzymes\t{describe_enzymes(enzymes)}\n")
        fh.write(f"concatemers\t{stats.n_concatemers}\n")
        fh.write(f"monomers\t{stats.n_monomers}\n")
        fh.write(f"cut_points\t{stats.n_cut_points}\n")
        fh.write(f"shared_cut_points\t{stats.n_shared_cut_points}\n")
        for enz in enzymes:
            n = stats.cuts_per_enzyme.get(enz.name, 0)
            fh.write(f"cuts_{enz.name}\t{n}\n")


def run_digest(args) -> int:
    try:
        enzymes = resolve_enzymes(args.enzyme)
    except EnzymeSpecError as exc:
        logger.error("%s", exc)
        return 2

    logger.info("Digesting with %s", describe_enzymes(enzymes))
    if args.dry_run:
        for enz in enzymes:
            print(f"{enz.name}\t{enz.site}\tfst5={enz.fst5}\t"
                  f"palindromic={enz.is_palindromic}")
        return 0

    stats = DigestStats()
    infile = _open_input(args.input, args.header)
    try:
        header = _output_header(infile, args.header, enzymes)
        mode = "wb" if args.output != "-" else "wb0"
        with pysam.AlignmentFile(
            args.output, mode, header=header,
            threads=max(1, args.threads),
        ) as out:
            reads = get_concatemer_seqs(
                infile, enzymes, remove_tags=args.remove_tags, stats=stats)
            for read in reads:
                if args.max_reads and stats.n_concatemers > args.max_reads:
                    break
                out.write(read)
    finally:
        infile.close()

    for line in stats.summary_lines(enzymes):
        logger.info("%s", line)
    if args.stats:
        _write_stats(args.stats, enzymes, stats)
        logger.info("Wrote per-enzyme report to %s", args.stats)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    _setup_logging(getattr(args, "log_level", logging.INFO),
                   getattr(args, "logfile", None))
    if args.command == "digest":
        return run_digest(args)
    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
