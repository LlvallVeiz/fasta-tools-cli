"""Record operations: filter, extract, per-record GC, reverse complement,
and format conversion.
"""

from __future__ import annotations

import re

from .parser import FastaRecord, gc_content, gc_windows, reverse_complement


def filter_records(
    records: list[FastaRecord],
    min_length: int | None = None,
    max_length: int | None = None,
    id_pattern: str | None = None,
    header_pattern: str | None = None,
    use_regex: bool = False,
) -> list[FastaRecord]:
    """Keep records matching all given criteria.

    Length bounds are inclusive. Patterns are matched as case-insensitive
    substrings by default; with ``use_regex=True`` they are treated as
    regular expressions (also case-insensitive). ``id_pattern`` matches
    against the record ID only, ``header_pattern`` against the full
    definition line (ID + description).
    """
    def matches(pattern: str, text: str) -> bool:
        if use_regex:
            return re.search(pattern, text, re.IGNORECASE) is not None
        return pattern.lower() in text.lower()

    kept: list[FastaRecord] = []
    for record in records:
        if min_length is not None and len(record) < min_length:
            continue
        if max_length is not None and len(record) > max_length:
            continue
        if id_pattern is not None and not matches(id_pattern, record.record_id):
            continue
        if header_pattern is not None and not matches(header_pattern, record.header):
            continue
        kept.append(record)
    return kept


def extract_records(
    records: list[FastaRecord],
    wanted_ids: list[str],
) -> tuple[list[FastaRecord], list[str]]:
    """Return records whose ID is in *wanted_ids*, in requested order.

    Returns ``(found, missing)``: ``found`` follows the order of
    *wanted_ids* (duplicated IDs in the request are returned once);
    ``missing`` lists requested IDs that do not occur in *records*.
    """
    by_id = {r.record_id: r for r in records}
    found: list[FastaRecord] = []
    missing: list[str] = []
    for wanted in dict.fromkeys(wanted_ids):  # de-duplicate, keep order
        record = by_id.get(wanted)
        if record is None:
            missing.append(wanted)
        else:
            found.append(record)
    return found, missing


def per_record_gc(records: list[FastaRecord]) -> list[tuple[str, int, float]]:
    """Return ``(record_id, length, gc_percent)`` per record."""
    return [(r.record_id, len(r), gc_content(r.sequence)) for r in records]


def window_gc_table(
    records: list[FastaRecord],
    window: int,
    step: int | None = None,
) -> list[tuple[str, int, float]]:
    """Return ``(record_id, window_start_1based, gc_percent)`` rows.

    One row per window per record; starts are 1-based positions.
    """
    rows: list[tuple[str, int, float]] = []
    for record in records:
        for i, value in enumerate(
            gc_windows(record.sequence, window, step)
        ):
            start = i * (window if step is None else step) + 1
            rows.append((record.record_id, start, value))
    return rows


def revcomp_records(records: list[FastaRecord]) -> list[FastaRecord]:
    """Return new records with reverse-complemented sequences.

    IDs are kept; the description gains a ``[reverse-complement]``
    marker so output files stay traceable to their source record.
    """
    out: list[FastaRecord] = []
    for record in records:
        desc = (record.description + " " if record.description else "")
        desc += "[reverse-complement]"
        out.append(
            FastaRecord(
                record_id=record.record_id,
                description=desc,
                sequence=reverse_complement(record.sequence),
            )
        )
    return out


def to_tsv(records: list[FastaRecord]) -> str:
    """Summarise records as a TSV table: id, length, gc%, description."""
    lines = ["id\tlength\tgc_percent\tdescription"]
    for record in records:
        lines.append(
            f"{record.record_id}\t{len(record)}\t"
            f"{gc_content(record.sequence):.2f}\t{record.description}"
        )
    return "\n".join(lines) + "\n"


def to_oneline_fasta(records: list[FastaRecord]) -> str:
    """Serialise records as FASTA with one sequence line per record."""
    return "".join(r.to_fasta(line_width=10 ** 9) for r in records)
