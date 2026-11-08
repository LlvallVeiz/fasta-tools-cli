"""fasta-tools-cli: a multi-purpose FASTA toolkit in pure Python."""

__version__ = "1.0.0"

from .parser import FastaError, FastaRecord, parse_fasta, parse_fasta_text
from .stats import FastaStats, compute_n50_l50, summarise

__all__ = [
    "__version__",
    "FastaError",
    "FastaRecord",
    "FastaStats",
    "compute_n50_l50",
    "parse_fasta",
    "parse_fasta_text",
    "summarise",
]
