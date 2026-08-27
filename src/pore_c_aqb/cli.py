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
import contextlib
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
from pore_c_aqb.report import (
    describe_cut,
    describe_overhang,
    enzyme_table,
    expected_site_spacing,
)

PROG = "pore-c-aqb"
logger = logging.getLogger(PROG)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description=(
            "Pore-C / CiFi tools: a modified pore-c-py that digests with "
            "several enzymes at once, undoes the cuts the enzyme never made, "
            "and tells you which enzymes actually cut your library."
        ),
        epilog=(
            "  pore-c-aqb digest DpnII,NlaIII reads.bam --output mono.bam\n"
            "  pore-c-aqb merge aligned.ns.bam --output fragments.tsv.gz\n"
            "  pore-c-aqb junctions aligned.ns.bam ref.fa --enzymes "
            "DpnII,NlaIII\n"
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
        "--version", action="version", version=f"{PROG} {__version__}")
    d.add_argument(
        "enzyme",
        help=("One or more restriction enzyme names. Separate several with "
              "commas, e.g. 'DpnII,NlaIII'. Names are those of Biopython's "
              "Restriction module."),
    )
    d.add_argument(
        "input", nargs="*", default=None,
        help=("Unaligned BAM(s) of concatemers, or a directory, or '-' for "
              "stdin. wf-pore-c calls the digest with the file before the "
              "enzyme in one of its two branches, so that order is accepted "
              "too."))
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
        "--max_reads", type=int, default=None,
        help="Take only the first N concatemers. Useful for testing.")
    d.add_argument(
        "--max_monomers", type=int, default=None,
        help=("Drop a concatemer cut into more than this many monomers. "
              "The whole read is excluded, not trimmed."))
    d.add_argument(
        "--excluded_list", type=Path, default=None,
        help="Write the name of each excluded read to this file.")
    d.add_argument(
        "--excluded_bam", type=Path, default=None,
        help="Write each excluded read to this BAM.")
    d.add_argument(
        "--recursive", action="store_true",
        help="If INPUT is a directory, search it recursively.")
    d.add_argument(
        "--glob", default="*.bam",
        help="If INPUT is a directory, match files with this glob.")
    d.add_argument(
        "--threads", type=int, default=1,
        help="Compute threads for BAM compression.")
    d.add_argument(
        "--stats", type=Path, default=None,
        help="Write the per-enzyme site report to this file as TSV.")
    d.add_argument(
        "--dry-run", action="store_true",
        help=("Print what each enzyme will do, then exit without reading "
              "input. Use it to check the enzymes before a long run."))

    sub.add_parser(
        "merge", add_help=False,
        help=("Glue back fragments the in silico digest split but the enzyme "
              "never cut. Fixes the inflated Hi-C diagonal."))
    sub.add_parser(
        "junctions", add_help=False,
        help=("Which enzymes actually cut this library? Motif enrichment at "
              "ligation junctions, from aligned monomers."))

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


def split_positionals(enzyme: str, inputs: list[str]) -> tuple[str, list[str]]:
    """Work out which positional is the enzyme.

    wf-pore-c calls the digest both ways round: ``digest "$cutter"`` reading
    stdin in its chunked branch, and ``digest concatemers.bam "$cutter"`` in
    the other. Upstream's parser puts the enzyme last; this one puts it first
    so that ``digest DpnII,NlaIII reads.bam`` reads naturally. Accepting both
    costs one check and avoids a confusing failure.
    """
    inputs = list(inputs or [])
    if inputs and Path(enzyme).exists() and not Path(inputs[-1]).exists():
        # first positional is a real file, last one is not: enzyme is last
        return inputs[-1], [enzyme] + inputs[:-1]
    return enzyme, inputs


def resolve_inputs(inputs: list[str], glob: str, recursive: bool) -> list[str]:
    """Expand directories, and default to stdin when nothing was given."""
    if not inputs:
        return ["-"]
    found: list[str] = []
    for item in inputs:
        path = Path(item)
        if item != "-" and path.is_dir():
            matches = sorted(
                path.rglob(glob) if recursive else path.glob(glob))
            if not matches:
                raise SystemExit(
                    f"No file matching {glob!r} in directory {item}.")
            found.extend(str(m) for m in matches)
        else:
            found.append(item)
    return found


def _open_input(path: str, header_path: Path | None):
    if path == "-":
        if header_path is None:
            raise SystemExit(
                "Reading from stdin requires --header: a piped BAM carries no "
                "header that can be copied to the output."
            )
        return pysam.AlignmentFile("-", "rb", check_sq=False)
    if not Path(path).exists():
        raise SystemExit(f"Input file not found: {path}")
    try:
        return pysam.AlignmentFile(path, "rb", check_sq=False)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Cannot read {path} as a BAM: {exc}")


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
        fh.write(f"excluded_concatemers\t{stats.n_excluded}\n")
        fh.write(f"cut_points\t{stats.n_cut_points}\n")
        fh.write(f"shared_cut_points\t{stats.n_shared_cut_points}\n")
        fh.write(f"bases_read\t{stats.n_bases}\n")
        for enz in enzymes:
            n = stats.sites_per_enzyme.get(enz.name, 0)
            fh.write(f"site_{enz.name}\t{enz.site}\n")
            fh.write(f"cut_{enz.name}\t{describe_cut(enz)}\n")
            fh.write(f"overhang_{enz.name}\t{describe_overhang(enz)}\n")
            fh.write(f"sites_found_{enz.name}\t{n}\n")
            fh.write(f"pct_of_cut_points_{enz.name}\t"
                     f"{100.0 * n / stats.n_cut_points:.2f}\n"
                     if stats.n_cut_points else
                     f"pct_of_cut_points_{enz.name}\t0.00\n")
            obs = stats.n_bases / n if n else 0
            fh.write(f"observed_spacing_bp_{enz.name}\t{obs:.0f}\n")
            fh.write(f"chance_spacing_bp_{enz.name}\t"
                     f"{expected_site_spacing(enz.site):.0f}\n")


def run_digest(args) -> int:
    enzyme_spec, inputs = split_positionals(args.enzyme, args.input)
    try:
        enzymes = resolve_enzymes(enzyme_spec)
    except EnzymeSpecError as exc:
        logger.error("%s", exc)
        return 2

    if args.dry_run:
        # stdout, not the logger: --dry-run must still print under --quiet,
        # and its output is meant to be readable and pipeable.
        print(f"Enzymes resolved from {enzyme_spec!r}:")
        for line in enzyme_table(enzymes):
            print(line)
        print("Dry run: no data read. Remove --dry-run to digest.")
        return 0

    inputs = resolve_inputs(inputs, args.glob, args.recursive)
    logger.info("Digesting with %s", describe_enzymes(enzymes))
    for line in enzyme_table(enzymes):
        logger.info("%s", line)

    stats = DigestStats()
    mode = "wb" if args.output != "-" else "wb0"
    remaining = args.max_reads

    with contextlib.ExitStack() as stack:
        # stdin can only be opened once: read the header from this very handle
        # and keep it open, rather than reopening the stream later
        first = stack.enter_context(_open_input(inputs[0], args.header))
        header = _output_header(first, args.header, enzymes)

        out = stack.enter_context(pysam.AlignmentFile(
            args.output, mode, header=header, threads=max(1, args.threads)))
        excluded_list = (stack.enter_context(open(args.excluded_list, "w"))
                         if args.excluded_list else None)
        excluded_bam = (stack.enter_context(pysam.AlignmentFile(
            str(args.excluded_bam), "wb", header=header,
            threads=max(1, args.threads))) if args.excluded_bam else None)

        def on_excluded(read):
            if excluded_list is not None:
                excluded_list.write(f"{read.query_name}\n")
            if excluded_bam is not None:
                excluded_bam.write(read)

        for index, path in enumerate(inputs):
            if remaining is not None and remaining <= 0:
                break
            infile = first if index == 0 else stack.enter_context(
                _open_input(path, args.header))
            before = stats.n_concatemers
            for read in get_concatemer_seqs(
                    infile, enzymes, remove_tags=args.remove_tags, stats=stats,
                    max_monomers=args.max_monomers, on_excluded=on_excluded,
                    max_reads=remaining):
                out.write(read)
            if remaining is not None:
                remaining -= stats.n_concatemers - before

    for line in stats.summary_lines(enzymes):
        logger.info("%s", line)
    if args.stats:
        _write_stats(args.stats, enzymes, stats)
        logger.info("Wrote per-enzyme report to %s", args.stats)
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # merge and junctions own their (long) help text, so hand the rest of the
    # command line straight to them rather than re-declaring it here
    if argv and argv[0] == "merge":
        from pore_c_aqb.merge import main as merge_main
        return merge_main(argv[1:])
    if argv and argv[0] == "junctions":
        from pore_c_aqb.junctions import main as junctions_main
        return junctions_main(argv[1:])

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
