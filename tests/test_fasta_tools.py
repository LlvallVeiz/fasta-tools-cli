"""Tests for fasta-tools-cli (unittest, standard library only).

Every numeric expectation below is computed BY HAND in the comments —
no value is copied from program output. Run with:

    python -m unittest discover -s tests -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from fasta_tools import parse_fasta_text  # noqa: E402
from fasta_tools.ops import (  # noqa: E402
    extract_records,
    filter_records,
    per_record_gc,
    revcomp_records,
    to_oneline_fasta,
    to_tsv,
)
from fasta_tools.parser import (  # noqa: E402
    FastaError,
    FastaRecord,
    gc_content,
    gc_windows,
    reverse_complement,
)
from fasta_tools.stats import compute_n50_l50, summarise  # noqa: E402


class TestParser(unittest.TestCase):
    def test_multiline_lowercase_and_blank_lines(self):
        # "acgt" (lowercase, 4 bp) + "ACGT" (4 bp), blank line between:
        # joined and uppercased -> 8 bp "ACGTACGT".
        text = ">seq1 some description\nacgt\n\nACGT\n"
        records = parse_fasta_text(text)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].record_id, "seq1")
        self.assertEqual(records[0].description, "some description")
        self.assertEqual(records[0].sequence, "ACGTACGT")

    def test_header_without_description(self):
        records = parse_fasta_text(">only_id\nACGT\n")
        self.assertEqual(records[0].record_id, "only_id")
        self.assertEqual(records[0].description, "")

    def test_multiple_records_in_order(self):
        records = parse_fasta_text(">a\nAAA\n>b\nTTT\n>c\nGGG\n")
        self.assertEqual([r.record_id for r in records], ["a", "b", "c"])
        self.assertEqual([r.sequence for r in records], ["AAA", "TTT", "GGG"])

    def test_empty_file_raises(self):
        for text in ("", "   \n\n  \n"):
            with self.assertRaises(FastaError):
                parse_fasta_text(text)

    def test_sequence_before_header_raises(self):
        with self.assertRaises(FastaError):
            parse_fasta_text("ACGT\n>seq1\nACGT\n")

    def test_record_without_sequence_raises(self):
        with self.assertRaises(FastaError):
            parse_fasta_text(">seq1\n>seq2\nACGT\n")

    def test_invalid_character_raises(self):
        # 'J' and '1' are not IUPAC nucleotide codes.
        with self.assertRaises(FastaError):
            parse_fasta_text(">seq1\nACGTJ1\n")

    def test_duplicate_id_raises(self):
        with self.assertRaises(FastaError):
            parse_fasta_text(">dup\nACGT\n>dup\nTTTT\n")

    def test_empty_definition_line_raises(self):
        with self.assertRaises(FastaError):
            parse_fasta_text(">\nACGT\n")

    def test_ambiguous_n_accepted(self):
        records = parse_fasta_text(">n1\nACGTN\n")
        self.assertEqual(records[0].sequence, "ACGTN")

    def test_semicolon_comment_lines_skipped(self):
        records = parse_fasta_text(">s1\n; a Pearson FASTA comment\nACGT\n")
        self.assertEqual(records[0].sequence, "ACGT")


class TestGcContent(unittest.TestCase):
    def test_gc_percent_hand_computed(self):
        # "ACGT": G+C = 2, length 4 -> 50.0%
        self.assertEqual(gc_content("ACGT"), 50.0)
        # "GGCC": G+C = 4, length 4 -> 100.0%
        self.assertEqual(gc_content("GGCC"), 100.0)
        # "AATT": G+C = 0, length 4 -> 0.0%
        self.assertEqual(gc_content("AATT"), 0.0)
        # "ACGTN": G+C = 2, length 5 (N in denominator) -> 40.0%
        self.assertEqual(gc_content("ACGTN"), 40.0)
        # lowercase input is uppercased first: "ggcc" -> 100.0%
        self.assertEqual(gc_content("ggcc"), 100.0)
        # empty sequence -> 0.0 (defined, no division by zero)
        self.assertEqual(gc_content(""), 0.0)

    def test_gc_windows_hand_computed(self):
        # "GGGGAAAA", window 4, step 4: [GGGG -> 100%, AAAA -> 0%]
        self.assertEqual(gc_windows("GGGGAAAA", 4), [100.0, 0.0])
        # "GGCCAA", window 4, step 2: windows GGCC (100%), CCAA (50%),
        # AA (0%, trailing partial window of 2 bp)
        self.assertEqual(gc_windows("GGCCAA", 4, step=2), [100.0, 50.0, 0.0])

    def test_gc_windows_invalid_args(self):
        with self.assertRaises(ValueError):
            gc_windows("ACGT", 0)
        with self.assertRaises(ValueError):
            gc_windows("ACGT", 2, step=0)


class TestReverseComplement(unittest.TestCase):
    def test_revcomp_hand_computed(self):
        # "AAAACCCC": reverse -> "CCCCAAAA", complement -> "GGGGTTTT"
        self.assertEqual(reverse_complement("AAAACCCC"), "GGGGTTTT")
        # "ATGC": reverse -> "CGTA", complement -> "GCAT"
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_revcomp_ambiguous_bases(self):
        # "ARYN": reverse -> "NYRA"; complement N->N, Y->R, R->Y, A->T
        # -> "NRYT"
        self.assertEqual(reverse_complement("ARYN"), "NRYT")
        # Single-base check for every IUPAC code against its complement.
        pairs = {"A": "T", "T": "A", "C": "G", "G": "C",
                 "R": "Y", "Y": "R", "S": "S", "W": "W",
                 "K": "M", "M": "K", "B": "V", "V": "B",
                 "D": "H", "H": "D", "N": "N"}
        for base, comp in pairs.items():
            self.assertEqual(reverse_complement(base), comp)


class TestN50L50(unittest.TestCase):
    def test_n50_l50_hand_computed(self):
        # Lengths [10, 20, 30, 40]: total 100, half = 50.
        # Sorted desc: [40, 30, 20, 10]; cumulative: 40, 70, 90, 100.
        # First cumulative >= 50 is 70, at the 2nd record (30 bp).
        # -> N50 = 30, L50 = 2.
        self.assertEqual(compute_n50_l50([10, 20, 30, 40]), (30, 2))
        # Single record: N50 is its length, L50 = 1.
        self.assertEqual(compute_n50_l50([42]), (42, 1))
        # Empty list -> (0, 0).
        self.assertEqual(compute_n50_l50([]), (0, 0))
        # Equal lengths [25, 25, 25, 25]: total 100, half 50;
        # cumulative 25, 50 -> reaches half at 2nd record -> (25, 2).
        self.assertEqual(compute_n50_l50([25, 25, 25, 25]), (25, 2))


class TestSummarise(unittest.TestCase):
    def test_summarise_hand_computed(self):
        # Records: "ACGT" (4 bp, GC 50%), "GGGGCCCC" (8 bp, GC 100%),
        # "AA" (2 bp, GC 0%).
        # total = 14, min = 2, max = 8, mean = 14/3 = 4.6667,
        # median = 4, GC overall = (2+8+0)/14 = 10/14 = 71.4286%,
        # lengths desc [8, 4, 2], half = 7, cum 8 >= 7 at 1st -> N50=8, L50=1.
        records = [
            FastaRecord("r1", "", "ACGT"),
            FastaRecord("r2", "", "GGGGCCCC"),
            FastaRecord("r3", "", "AA"),
        ]
        stats = summarise(records)
        self.assertEqual(stats.n_records, 3)
        self.assertEqual(stats.total_length, 14)
        self.assertEqual(stats.min_length, 2)
        self.assertEqual(stats.max_length, 8)
        self.assertAlmostEqual(stats.mean_length, 14 / 3, places=4)
        self.assertEqual(stats.median_length, 4.0)
        self.assertAlmostEqual(stats.gc_percent, 100 * 10 / 14, places=4)
        self.assertEqual(stats.n50, 8)
        self.assertEqual(stats.l50, 1)
        self.assertEqual(stats.n_bases, 0)
        self.assertEqual(stats.ambiguous_bases, 0)

    def test_summarise_counts_n_and_ambiguous(self):
        # "ACGTN" (1 N) + "ACGR" (1 R, ambiguous non-N)
        records = [
            FastaRecord("a", "", "ACGTN"),
            FastaRecord("b", "", "ACGR"),
        ]
        stats = summarise(records)
        self.assertEqual(stats.n_bases, 1)
        self.assertEqual(stats.ambiguous_bases, 1)

    def test_summarise_empty_raises(self):
        with self.assertRaises(ValueError):
            summarise([])


class TestOps(unittest.TestCase):
    def setUp(self):
        self.records = [
            FastaRecord("short", "a short one", "ACGT"),          # 4 bp
            FastaRecord("medium_gc", "GC rich", "GGGGCCCC"),      # 8 bp
            FastaRecord("long_at", "AT rich", "AAAATTTTAAAA"),    # 12 bp
        ]

    def test_filter_by_length(self):
        kept = filter_records(self.records, min_length=5)
        self.assertEqual([r.record_id for r in kept],
                         ["medium_gc", "long_at"])
        kept = filter_records(self.records, max_length=8)
        self.assertEqual([r.record_id for r in kept],
                         ["short", "medium_gc"])
        kept = filter_records(self.records, min_length=5, max_length=8)
        self.assertEqual([r.record_id for r in kept], ["medium_gc"])

    def test_filter_by_pattern(self):
        kept = filter_records(self.records, id_pattern="GC")
        self.assertEqual([r.record_id for r in kept], ["medium_gc"])
        # case-insensitive substring on the full header
        kept = filter_records(self.records, header_pattern="at rich")
        self.assertEqual([r.record_id for r in kept], ["long_at"])
        # regex mode: IDs ending in "_at"
        kept = filter_records(self.records, id_pattern=r"_at$",
                              use_regex=True)
        self.assertEqual([r.record_id for r in kept], ["long_at"])

    def test_extract_order_and_missing(self):
        found, missing = extract_records(
            self.records, ["long_at", "short", "nope"])
        # found follows the REQUESTED order, not file order
        self.assertEqual([r.record_id for r in found],
                         ["long_at", "short"])
        self.assertEqual(missing, ["nope"])

    def test_per_record_gc(self):
        # short "ACGT" -> 50%, medium "GGGGCCCC" -> 100%,
        # long "AAAATTTTAAAA" -> 0%
        rows = per_record_gc(self.records)
        self.assertEqual(rows, [
            ("short", 4, 50.0),
            ("medium_gc", 8, 100.0),
            ("long_at", 12, 0.0),
        ])

    def test_revcomp_records(self):
        out = revcomp_records([self.records[0]])  # "ACGT" -> "ACGT"
        self.assertEqual(out[0].sequence, "ACGT")
        self.assertIn("[reverse-complement]", out[0].description)
        self.assertEqual(out[0].record_id, "short")
        out = revcomp_records([self.records[1]])  # "GGGGCCCC" -> "GGGGCCCC"
        self.assertEqual(out[0].sequence, "GGGGCCCC")
        out = revcomp_records([FastaRecord("x", "", "AAAACCCC")])
        self.assertEqual(out[0].sequence, "GGGGTTTT")

    def test_to_tsv(self):
        tsv = to_tsv(self.records)
        lines = tsv.strip().split("\n")
        self.assertEqual(lines[0], "id\tlength\tgc_percent\tdescription")
        self.assertEqual(lines[1], "short\t4\t50.00\ta short one")
        self.assertEqual(lines[2], "medium_gc\t8\t100.00\tGC rich")

    def test_to_oneline_fasta(self):
        text = to_oneline_fasta([self.records[0]])
        self.assertEqual(text, ">short a short one\nACGT\n")


class TestSampleFile(unittest.TestCase):
    """End-to-end check against examples/sample.fasta (real file)."""

    def test_sample_stats_hand_computed(self):
        # sample.fasta: 5 records, lengths 20/30/40/15/60 = 165 bp total.
        # GC counts: 10 + 30 + 0 + 6 + 17 = 63 -> 63/165 = 38.1818%.
        # N50: sorted desc [60, 40, 30, 20, 15], half = 82.5;
        # cumulative 60, then 100 >= 82.5 at 2nd record -> N50 = 40, L50 = 2.
        # N bases: 3 (in syn_ambig "ACGTN" x3).
        path = os.path.join(os.path.dirname(__file__), "..",
                            "examples", "sample.fasta")
        from fasta_tools.parser import parse_fasta  # noqa: E402
        records = parse_fasta(path)
        stats = summarise(records)
        self.assertEqual(stats.n_records, 5)
        self.assertEqual(stats.total_length, 165)
        self.assertEqual(stats.min_length, 15)
        self.assertEqual(stats.max_length, 60)
        self.assertAlmostEqual(stats.mean_length, 33.0)
        self.assertEqual(stats.median_length, 30.0)
        self.assertAlmostEqual(stats.gc_percent, 100 * 63 / 165, places=2)
        self.assertEqual(stats.n50, 40)
        self.assertEqual(stats.l50, 2)
        self.assertEqual(stats.n_bases, 3)


class TestCli(unittest.TestCase):
    """Smoke tests for the CLI entry point (subprocess-free)."""

    def _run_cli(self, argv):
        from fasta_tools.cli import main  # noqa: E402
        return main(argv)

    def test_cli_stats_and_gc_exit_zero(self):
        path = os.path.join(os.path.dirname(__file__), "..",
                            "examples", "sample.fasta")
        self.assertEqual(self._run_cli(["stats", path]), 0)
        self.assertEqual(self._run_cli(["gc", path]), 0)
        self.assertEqual(self._run_cli(["revcomp", path]), 0)
        self.assertEqual(self._run_cli(["convert", "--to", "tsv", path]), 0)
        self.assertEqual(
            self._run_cli(["convert", "--to", "oneline", path]), 0)

    def test_cli_extract_missing_id_exits_one(self):
        path = os.path.join(os.path.dirname(__file__), "..",
                            "examples", "sample.fasta")
        self.assertEqual(
            self._run_cli(["extract", "--id", "does_not_exist", path]), 1)

    def test_cli_extract_from_id_file(self):
        base = os.path.dirname(__file__)
        fasta = os.path.join(base, "..", "examples", "sample.fasta")
        ids = os.path.join(base, "..", "examples", "ids.txt")
        with tempfile.NamedTemporaryFile("w", suffix=".fasta",
                                        delete=False) as tmp:
            out_path = tmp.name
        try:
            self.assertEqual(
                self._run_cli(["extract", "--id-file", ids,
                               "-o", out_path, fasta]), 0)
            with open(out_path, encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn(">syn_gc_rich", text)
            self.assertIn(">syn_ambig", text)
            self.assertNotIn(">syn_short", text)
        finally:
            os.unlink(out_path)


if __name__ == "__main__":
    unittest.main()
