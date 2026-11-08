# fasta-tools-cli

A multi-purpose FASTA toolkit in pure Python — no external dependencies.
Parse, summarise, filter, extract, and transform FASTA files from one
small command-line tool.

```text
$ ./fasta-tools stats examples/sample.fasta
records:        5
total length:   165 bp
min length:     15 bp
max length:     60 bp
mean length:    33.00 bp
median length:  30.00 bp
GC content:     38.18%
N50:            40 bp
L50:            2
N bases:        3
ambiguous:      0
```

## Background

FASTA is the plain-text format bioinformatics uses to store nucleotide
and protein sequences: each record is a one-line *definition line*
starting with `>` (an identifier plus a free-text description),
followed by one or more lines of sequence data. The format is specified
by NCBI; sequences use the IUPAC single-letter nucleotide codes,
including `N` for an unknown base and the other ambiguity codes
(`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`).

Two derived quantities deserve a definition, because tools disagree on
them unless stated:

- **GC%** = `100 * (G + C) / sequence length`. Ambiguous bases count
  in the denominator but not in the numerator — GC% is reported over
  raw sequence length.
- **N50 / L50**: sort the record lengths descending and accumulate;
  N50 is the length at which the running sum first reaches half of the
  total length, and L50 is how many records that took. (For the example
  file: lengths 60, 40, 30, 20, 15 — half of 165 is 82.5; 60 then
  60+40 = 100 ≥ 82.5, so N50 = 40 and L50 = 2.)

The full research notes — FASTA spec, the IUPAC table with complements,
N50/L50, and prior art (seqtk, Biopython, seqkit) — are in
[`docs/RESEARCH.md`](docs/RESEARCH.md).

## Installation

Requires Python 3.9+ and nothing else (standard library only).

```bash
git clone https://github.com/LlvallVeiz/fasta-tools-cli.git
cd fasta-tools-cli
```

No `pip install` step: run `./fasta-tools ...` directly, or
`PYTHONPATH=src python -m fasta_tools ...`.

## Usage

```
fasta-tools {stats, filter, extract, gc, revcomp, convert} [options] <input.fasta>
```

| Subcommand | What it does |
|---|---|
| `stats` | Record count, total/min/max/mean/median length, GC%, N50, L50, N-base and ambiguity counts |
| `filter` | Keep records by `--min-length` / `--max-length` and/or `--id-pattern` / `--header-pattern` (substring, or regex with `--regex`) |
| `extract` | Pull out records by ID: repeatable `--id` and/or `--id-file` (one ID per line) |
| `gc` | GC% per record; with `--window N` (and optional `--step`), sliding-window GC% |
| `revcomp` | Reverse complement every record (full IUPAC table) |
| `convert` | `--to tsv` (id / length / GC% / description) or `--to oneline` (one sequence line per record) |

Output goes to stdout unless `-o/--output FILE` is given (supported by
`filter`, `extract`, `revcomp`, `convert`). Errors are reported on
stderr with a non-zero exit code: exit 2 for unusable input (missing
file, empty file, malformed FASTA), exit 1 when `extract` is asked for
an ID that does not exist.

## Examples

All output below is real output from running the tool on
`examples/sample.fasta` (5 synthetic records, 165 bp total; see
[`examples/example_output.txt`](examples/example_output.txt) for the
full transcript).

```bash
# Summary statistics
./fasta-tools stats examples/sample.fasta

# Records of at least 30 bp
./fasta-tools filter --min-length 30 examples/sample.fasta
# -> syn_gc_rich (30 bp), syn_at_rich (40 bp), syn_mixed (60 bp)

# Records whose header mentions GC-rich
./fasta-tools filter --header-pattern GC-rich examples/sample.fasta
# -> syn_gc_rich only

# Extract two records listed in a file
./fasta-tools extract --id-file examples/ids.txt examples/sample.fasta
# -> syn_gc_rich, syn_ambig

# GC% per record
./fasta-tools gc examples/sample.fasta
# syn_short   20  50.00
# syn_gc_rich 30  100.00
# syn_at_rich 40  0.00
# syn_ambig   15  40.00
# syn_mixed   60  28.33

# Sliding-window GC% (10 bp windows)
./fasta-tools gc --window 10 examples/sample.fasta

# Reverse complement (note syn_ambig: ACGTNACGTNACGTN -> NACGTNACGTNACGT)
./fasta-tools revcomp examples/sample.fasta

# Convert to a TSV summary
./fasta-tools convert --to tsv examples/sample.fasta

# Convert to one-line-per-record FASTA
./fasta-tools convert --to oneline examples/sample.fasta
```

Error handling:

```bash
$ ./fasta-tools stats /dev/null
fasta-tools: error: /dev/null: file is empty (no FASTA records found)
# exit code 2

$ ./fasta-tools extract --id no_such_id examples/sample.fasta
fasta-tools: error: ID(s) not found: no_such_id
# exit code 1
```

## Project structure

```text
fasta-tools-cli/
├── fasta-tools                 # launcher: ./fasta-tools <subcommand> ...
├── src/fasta_tools/
│   ├── __init__.py             # public API + version
│   ├── __main__.py             # `python -m fasta_tools`
│   ├── parser.py               # FASTA parser, FastaRecord, GC, reverse complement
│   ├── stats.py                # summary statistics incl. N50/L50
│   ├── ops.py                  # filter / extract / GC tables / conversion
│   └── cli.py                  # argparse CLI, six subcommands
├── tests/
│   └── test_fasta_tools.py     # 31 unittest tests, hand-computed expectations
├── examples/
│   ├── sample.fasta            # 5 synthetic records (motif repeats, no real organism)
│   ├── ids.txt                 # two IDs, for the extract example
│   └── example_output.txt      # real output of every subcommand on sample.fasta
├── docs/
│   └── RESEARCH.md             # FASTA spec, IUPAC table, N50/L50, prior art, references
├── requirements.txt            # (empty — standard library only)
├── .gitignore
└── LICENSE                     # MIT
```

## Tests

```bash
python -m unittest discover -s tests -v
# Ran 31 tests — OK
```

Every numeric expectation in the tests is computed **by hand** and the
arithmetic is written out in comments next to the assertion — e.g.
N50/L50 for lengths `[10, 20, 30, 40]` (total 100, half 50; cumulative
40, 70 → N50 = 30, L50 = 2) and `revcomp("ARYN") = "NRYT"`. Nothing is
copied from program output. The suite also runs each CLI subcommand on
the real example file and checks the sample file's statistics
end-to-end (5 records, 165 bp, GC 38.18%, N50 = 40, L50 = 2).

## Design decisions

- **Pure Python, zero dependencies.** FASTA parsing is string work;
  pulling in Biopython would buy format breadth this tool does not
  need, and would hide the parsing logic this project exists to show.
- **Strict where ambiguity corrupts results, tolerant where the format
  is.** Blank lines, lowercase bases, multi-line sequences and `;`
  comment lines are accepted; an empty file, data before the first
  header, a record without sequence, an invalid character, or a
  duplicated ID is a hard error with the file and line number in the
  message.
- **Uppercase normalisation at parse time.** All sequences are stored
  uppercased, so every later operation (GC, revcomp, filter) sees one
  canonical form.
- **Full IUPAC complement table.** Reverse complement handles ambiguous
  bases (`R`↔`Y`, `K`↔`M`, `N`→`N`, …), not just A/C/G/T.
- **GC% over raw length.** Ambiguous bases sit in the denominator;
  the choice is documented in `parser.py`, `docs/RESEARCH.md`, and
  above, so it can be challenged.
- **`extract` preserves the requested order**, not file order, and
  reports every missing ID instead of failing on the first one.
- **Output is plain text on stdout** (FASTA stays FASTA, tables are
  TSV), so every subcommand composes with Unix pipes and with
  `seqtk`/`seqkit` workflows.

## Limitations & future work

- **In-memory parsing.** The whole file is read at once; a multi-GB
  FASTA will use roughly its own size in RAM. A streaming parser
  (generator-based, one record at a time) is the natural next step.
- **No FASTQ support.** Quality scores are out of scope; `seqtk` or
  Biopython cover that.
- **No sequence validation beyond the IUPAC alphabet.** `U` is accepted
  (RNA input) but the tool does not distinguish DNA from RNA, and it
  does not check that a sequence is a plausible coding sequence.
- **No gzip input.** Compressed `.fasta.gz` files must be decompressed
  first (`zcat file.fasta.gz | ...` is not supported either — input is
  a file path, not stdin).
- **No parallel processing.** All operations are single-threaded; fine
  for the file sizes this tool targets.
- **Substring matching is case-insensitive** in `filter`; exact-case
  matching is not offered.

## License

MIT — see [LICENSE](LICENSE).
