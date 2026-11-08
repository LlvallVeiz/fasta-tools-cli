"""Summary statistics over a set of FASTA records.

N50/L50 follow the standard assembly definition (see docs/RESEARCH.md):
sort the record lengths descending, walk the cumulative sum, and take
the length at which the cumulative sum first reaches at least half of
the total length. N50 is that length; L50 is how many records it took
to reach it.
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import mean, median

from .parser import FastaRecord, gc_content


@dataclass(frozen=True)
class FastaStats:
    """All summary statistics for one FASTA file."""

    n_records: int
    total_length: int
    min_length: int
    max_length: int
    mean_length: float
    median_length: float
    gc_percent: float
    n50: int
    l50: int
    n_bases: int  # count of 'N' (unknown) bases
    # Count of IUPAC codes other than A/C/G/T/U/N (e.g. R, Y, S, W...).
    # 'N' is reported separately on the n_bases line above.
    ambiguous_bases: int

    def as_lines(self) -> list[str]:
        """Human-readable ``key: value`` lines, in stable order."""
        return [
            f"records:        {self.n_records}",
            f"total length:   {self.total_length} bp",
            f"min length:     {self.min_length} bp",
            f"max length:     {self.max_length} bp",
            f"mean length:    {self.mean_length:.2f} bp",
            f"median length:  {self.median_length:.2f} bp",
            f"GC content:     {self.gc_percent:.2f}%",
            f"N50:            {self.n50} bp",
            f"L50:            {self.l50}",
            f"N bases:        {self.n_bases}",
            f"ambiguous:      {self.ambiguous_bases}",
        ]


def compute_n50_l50(lengths: list[int]) -> tuple[int, int]:
    """Return ``(N50, L50)`` for a list of sequence lengths.

    Hand-check example (also in the test suite): lengths
    ``[10, 20, 30, 40]`` total 100 bp, half = 50. Sorted descending
    ``[40, 30, 20, 10]``; cumulative sums are 40, 70, 90, 100, so
    the sum first reaches >= 50 at the second record (30 bp).
    Hence N50 = 30 and L50 = 2.
    """
    if not lengths:
        return 0, 0
    total = sum(lengths)
    half = total / 2
    running = 0
    for i, length in enumerate(sorted(lengths, reverse=True), start=1):
        running += length
        if running >= half:
            return length, i
    # Unreachable: the last cumulative sum equals the total >= half.
    return sorted(lengths)[0], len(lengths)


def summarise(records: list[FastaRecord]) -> FastaStats:
    """Compute summary statistics for a list of records.

    Raises:
        ValueError: if *records* is empty.
    """
    if not records:
        raise ValueError("cannot summarise an empty record list")
    lengths = [len(r) for r in records]
    total = sum(lengths)
    all_seq = "".join(r.sequence for r in records)
    n50, l50 = compute_n50_l50(lengths)
    return FastaStats(
        n_records=len(records),
        total_length=total,
        min_length=min(lengths),
        max_length=max(lengths),
        mean_length=mean(lengths),
        median_length=float(median(lengths)),
        gc_percent=gc_content(all_seq),
        n50=n50,
        l50=l50,
        n_bases=all_seq.count("N"),
        ambiguous_bases=sum(
            1 for base in all_seq if base not in "ACGTUN"
        ),
    )
