"""Command-line interface for fasta-tools-cli.

Six subcommands, all reading a FASTA file and writing to stdout unless
``-o/--output`` is given:

* ``stats``    — record count, length stats, GC%, N50, L50.
* ``filter``   — keep records by length and/or header/ID pattern.
* ``extract``  — pull out records by ID (``--id`` and/or ``--id-file``).
* ``gc``       — GC% per record, optionally in sliding windows.
* ``revcomp``  — reverse complement every record.
* ``convert``  — FASTA -> TSV summary or FASTA with one line per record.
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .ops import (
    extract_records,
    filter_records,
    per_record_gc,
    revcomp_records,
    to_oneline_fasta,
    to_tsv,
    window_gc_table,
)
from .parser import FastaError, parse_fasta
from .stats import summarise


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="fasta-tools",
        description=(
            "Multi-purpose FASTA toolkit: stats, filter, extract, GC, "
            "reverse complement, and format conversion. Pure Python, "
            "no external dependencies."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_input(p: argparse.ArgumentParser) -> None:
        p.add_argument("input", help="input FASTA file")

    def add_output(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "-o", "--output",
            help="write output to FILE instead of stdout",
        )

    # stats -----------------------------------------------------------------
    p_stats = sub.add_parser(
        "stats", help="summary statistics: counts, lengths, GC%%, N50, L50"
    )
    add_input(p_stats)
    p_stats.set_defaults(func=cmd_stats)

    # filter ----------------------------------------------------------------
    p_filter = sub.add_parser(
        "filter", help="keep records by length and/or ID/header pattern"
    )
    add_input(p_filter)
    add_output(p_filter)
    p_filter.add_argument("--min-length", type=int, default=None,
                          help="keep records with length >= N")
    p_filter.add_argument("--max-length", type=int, default=None,
                          help="keep records with length <= N")
    p_filter.add_argument("--id-pattern", default=None,
                          help="substring (or regex with --regex) in record ID")
    p_filter.add_argument("--header-pattern", default=None,
                          help="substring (or regex with --regex) in full header")
    p_filter.add_argument("--regex", action="store_true",
                          help="treat patterns as regular expressions")
    p_filter.set_defaults(func=cmd_filter)

    # extract ---------------------------------------------------------------
    p_extract = sub.add_parser(
        "extract", help="extract records by ID"
    )
    add_input(p_extract)
    add_output(p_extract)
    p_extract.add_argument("--id", action="append", default=[],
                           help="record ID to extract (repeatable)")
    p_extract.add_argument("--id-file", default=None,
                           help="file with one record ID per line")
    p_extract.set_defaults(func=cmd_extract)

    # gc --------------------------------------------------------------------
    p_gc = sub.add_parser(
        "gc", help="GC%% per record, optionally in sliding windows"
    )
    add_input(p_gc)
    p_gc.add_argument("--window", type=int, default=None,
                      help="window size in bp for sliding-window GC%%")
    p_gc.add_argument("--step", type=int, default=None,
                      help="step between windows (default: window size)")
    p_gc.set_defaults(func=cmd_gc)

    # revcomp ----------------------------------------------------------------
    p_revcomp = sub.add_parser(
        "revcomp", help="reverse complement every record"
    )
    add_input(p_revcomp)
    add_output(p_revcomp)
    p_revcomp.set_defaults(func=cmd_revcomp)

    # convert ----------------------------------------------------------------
    p_convert = sub.add_parser(
        "convert", help="convert FASTA to TSV summary or one-line FASTA"
    )
    add_input(p_convert)
    add_output(p_convert)
    p_convert.add_argument(
        "--to", choices=["tsv", "oneline"], required=True,
        help="target format: 'tsv' (id/length/GC%%/description) or "
             "'oneline' (FASTA, one sequence line per record)",
    )
    p_convert.set_defaults(func=cmd_convert)

    return parser


def _load(path: str):
    """Parse *path*, printing a clean error and exiting 2 on failure."""
    try:
        return parse_fasta(path)
    except FileNotFoundError:
        print(f"fasta-tools: error: file not found: {path}", file=sys.stderr)
        sys.exit(2)
    except FastaError as exc:
        print(f"fasta-tools: error: {exc}", file=sys.stderr)
        sys.exit(2)


def _emit(text: str, output: str | None) -> None:
    """Write *text* to *output* file or stdout."""
    if output:
        with open(output, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)


def _emit_fasta(records, output: str | None) -> None:
    _emit("".join(r.to_fasta() for r in records), output)


def cmd_stats(args: argparse.Namespace) -> int:
    records = _load(args.input)
    stats = summarise(records)
    print("\n".join(stats.as_lines()))
    return 0


def cmd_filter(args: argparse.Namespace) -> int:
    records = _load(args.input)
    kept = filter_records(
        records,
        min_length=args.min_length,
        max_length=args.max_length,
        id_pattern=args.id_pattern,
        header_pattern=args.header_pattern,
        use_regex=args.regex,
    )
    _emit_fasta(kept, args.output)
    print(f"fasta-tools: kept {len(kept)} of {len(records)} records",
          file=sys.stderr)
    return 0


def cmd_extract(args: argparse.Namespace) -> int:
    records = _load(args.input)
    wanted: list[str] = list(args.id)
    if args.id_file:
        try:
            with open(args.id_file, encoding="utf-8") as handle:
                wanted += [line.strip() for line in handle if line.strip()]
        except FileNotFoundError:
            print(f"fasta-tools: error: file not found: {args.id_file}",
                  file=sys.stderr)
            return 2
    if not wanted:
        print("fasta-tools: error: no IDs given "
              "(use --id and/or --id-file)", file=sys.stderr)
        return 2
    found, missing = extract_records(records, wanted)
    if missing:
        print("fasta-tools: error: ID(s) not found: " + ", ".join(missing),
              file=sys.stderr)
        return 1
    _emit_fasta(found, args.output)
    return 0


def cmd_gc(args: argparse.Namespace) -> int:
    records = _load(args.input)
    if args.window is not None:
        if args.window < 1 or (args.step is not None and args.step < 1):
            print("fasta-tools: error: --window and --step must be >= 1",
                  file=sys.stderr)
            return 2
        rows = window_gc_table(records, args.window, args.step)
        print("id\twindow_start\tgc_percent")
        for record_id, start, value in rows:
            print(f"{record_id}\t{start}\t{value:.2f}")
    else:
        print("id\tlength\tgc_percent")
        for record_id, length, value in per_record_gc(records):
            print(f"{record_id}\t{length}\t{value:.2f}")
    return 0


def cmd_revcomp(args: argparse.Namespace) -> int:
    records = _load(args.input)
    _emit_fasta(revcomp_records(records), args.output)
    return 0


def cmd_convert(args: argparse.Namespace) -> int:
    records = _load(args.input)
    if args.to == "tsv":
        _emit(to_tsv(records), args.output)
    else:
        _emit(to_oneline_fasta(records), args.output)
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point; returns a process exit code."""
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
