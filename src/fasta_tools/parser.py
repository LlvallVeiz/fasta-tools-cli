"""FASTA parser for fasta-tools-cli.

Implements the FASTA text format as described by NCBI: a record starts with
a single-line definition line beginning with ``>``, followed by one or more
lines of sequence data (NCBI recommends <= 80 characters per line).
The sequence identifier (SeqID) is the first whitespace-delimited token
after ``>``; the rest of the line is the description.

This parser is deliberately liberal where the format is liberal in
practice, and strict where ambiguity would silently corrupt results:

* multi-line sequences are joined;
* blank lines are skipped;
* lowercase bases are normalised to uppercase (as in the NCBI spec,
  lowercase letters are accepted and mapped to uppercase);
* ``;`` comment lines (original Pearson FASTA) are skipped;
* IUPAC ambiguity codes (including ``N``) are accepted;
* an empty file, a file with no records, sequence data before the first
  header, a record with no sequence, an invalid character, or a
  duplicated record ID all raise :class:`FastaError` with a clear message.
"""

from __future__ import annotations

from dataclasses import dataclass

# IUPAC nucleotide codes (Cornish-Bowden 1985; NCBI FASTA spec requires
# IUPAC symbols). Includes U for RNA input, which is kept as-is.
IUPAC_BASES = frozenset("ACGTUNRYSWKMBDHV")

# Complement table for reverse complement, including every IUPAC
# ambiguity code. Source: IUPAC-IUB nomenclature, Cornish-Bowden (1985),
# Nucl. Acids Res. 13: 3021-3030.
COMPLEMENT = str.maketrans({
    "A": "T", "T": "A", "U": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
})


class FastaError(ValueError):
    """Raised when a FASTA file cannot be parsed unambiguously."""


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record.

    Attributes:
        record_id: first whitespace-delimited token after ``>``.
        description: remainder of the definition line (may be empty).
        sequence: sequence string, uppercased, whitespace removed.
    """

    record_id: str
    description: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        """Full definition line without the leading ``>``."""
        if self.description:
            return f"{self.record_id} {self.description}"
        return self.record_id

    def to_fasta(self, line_width: int = 80) -> str:
        """Serialise the record as FASTA text (default 80 chars/line)."""
        lines = [f">{self.header}"]
        seq = self.sequence
        for i in range(0, len(seq), line_width):
            lines.append(seq[i:i + line_width])
        return "\n".join(lines) + "\n"


def parse_fasta_text(text: str, source: str = "<input>") -> list[FastaRecord]:
    """Parse FASTA records from a string.

    Args:
        text: full FASTA file contents.
        source: name shown in error messages (file path or ``<stdin>``).

    Returns:
        List of :class:`FastaRecord` in file order.

    Raises:
        FastaError: empty input, no records, data before the first
            header, record without sequence, invalid character, or
            duplicated record ID.
    """
    if not text.strip():
        raise FastaError(f"{source}: file is empty (no FASTA records found)")

    records: list[FastaRecord] = []
    seen_ids: set[str] = set()
    current_id: str | None = None
    current_desc = ""
    current_chunks: list[str] = []
    current_line = 0

    def flush(line_no: int) -> None:
        """Close the record currently being built, validating it."""
        nonlocal current_id, current_desc, current_chunks
        if current_id is None:
            return
        sequence = "".join(current_chunks)
        if not sequence:
            raise FastaError(
                f"{source}: record '{current_id}' (line {line_no}) "
                "has no sequence data"
            )
        if current_id in seen_ids:
            raise FastaError(
                f"{source}: duplicated record ID '{current_id}' "
                "(IDs must be unique)"
            )
        seen_ids.add(current_id)
        records.append(FastaRecord(current_id, current_desc, sequence))
        current_id = None
        current_desc = ""
        current_chunks = []

    for raw_line in text.splitlines():
        current_line += 1
        line = raw_line.strip()
        if not line:
            continue  # blank lines are ignored
        if line.startswith(";"):
            continue  # Pearson FASTA comment line
        if line.startswith(">"):
            flush(current_line)
            header = line[1:].strip()
            if not header:
                raise FastaError(
                    f"{source}: line {current_line}: empty definition line "
                    "('>' with no identifier)"
                )
            parts = header.split(None, 1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
            current_chunks = []
            continue
        if current_id is None:
            raise FastaError(
                f"{source}: line {current_line}: sequence data found before "
                "the first '>' definition line"
            )
        seq = line.upper().replace(" ", "").replace("\t", "")
        invalid = sorted(set(seq) - IUPAC_BASES)
        if invalid:
            raise FastaError(
                f"{source}: line {current_line}: invalid character(s) "
                f"{invalid} in sequence of record '{current_id}' "
                "(only IUPAC nucleotide codes are allowed)"
            )
        current_chunks.append(seq)

    flush(current_line)
    if not records:
        raise FastaError(f"{source}: no FASTA records found")
    return records


def parse_fasta(path: str) -> list[FastaRecord]:
    """Parse all records from a FASTA file at *path*.

    Raises:
        FastaError: if the file cannot be parsed (see
            :func:`parse_fasta_text`).
        FileNotFoundError: if *path* does not exist.
    """
    with open(path, "r", encoding="utf-8") as handle:
        return parse_fasta_text(handle.read(), source=path)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA/RNA sequence.

    Uses the full IUPAC complement table, so ambiguous bases are
    complemented correctly (e.g. ``R`` <-> ``Y``, ``N`` -> ``N``).
    The input is uppercased first; ``U`` complements to ``A``.
    """
    return sequence.upper().translate(COMPLEMENT)[::-1]


def gc_content(sequence: str) -> float:
    """GC percentage of *sequence*.

    Defined as ``100 * (G + C) / len(sequence)``; ambiguous bases
    (``N`` and other IUPAC codes) count in the denominator but not in
    the numerator. This is the convention used by most QC tools when
    reporting GC% over raw sequence length. Returns 0.0 for an empty
    sequence.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def gc_windows(sequence: str, window: int, step: int | None = None) -> list[float]:
    """Sliding-window GC% along *sequence*.

    Args:
        sequence: sequence to scan.
        window: window size in bases (must be >= 1).
        step: step between window starts; defaults to *window*
            (non-overlapping windows).

    Returns:
        GC% of each window. A trailing partial window shorter than
        *window* is included only if it is at least 1 base long and
        is not a duplicate of the previous window start; windows are
        taken while the start index is < len(sequence).

    Raises:
        ValueError: if *window* < 1 or *step* < 1.
    """
    if window < 1:
        raise ValueError("window size must be >= 1")
    step = window if step is None else step
    if step < 1:
        raise ValueError("step must be >= 1")
    seq = sequence.upper()
    return [gc_content(seq[i:i + window]) for i in range(0, len(seq), step)]
