# Research notes — fasta-tools-cli

Background research for the design decisions in this project. Every
claim below is traceable to a source that was actually read on
2026-11-08; nothing here is invented.

## 1. The FASTA format (NCBI)

The FASTA format is a text format for nucleotide or protein sequences:
a record starts with a single-line *definition line* beginning with
`>`, followed by raw sequence data. The SeqID is the first token after
`>` and must be unique per record; the rest of the line is a free-text
description. Sequence lines are recommended to be no longer than 80
characters, and only IUPAC symbols may be used in the sequence.

Design consequences for the parser in `src/fasta_tools/parser.py`:

- The record ID is the first whitespace-delimited token; the rest of
  the definition line is kept as `description`.
- Sequence lines are joined (multi-line sequences are normal).
- Lowercase letters are accepted and mapped to uppercase, as the NCBI
  specification allows.
- Blank lines are skipped (tolerated; many real-world files contain
  them).
- Lines starting with `;` are skipped: the original Pearson FASTA
  format allowed `;` comment lines, although NCBI tools do not
  recognise them. Skipping them is the tolerant choice.
- A file with no records, sequence data before the first `>`, a record
  with no sequence, an invalid character, or a duplicated ID is a hard
  error — silently guessing in these cases would corrupt downstream
  results.

Sources:

- NCBI, "FASTA Format for Nucleotide Sequences",
  https://www.ncbi.nlm.nih.gov/genbank/fastaformat/
- NCBI Bookshelf, GenBank Submissions Handbook, "FASTA Formatting",
  https://www.ncbi.nlm.nih.gov/books/NBK566986/
- Wikipedia, "FASTA format",
  https://en.wikipedia.org/wiki/FASTA_format (overview of the
  original Pearson format, including the `;` comment convention)

## 2. IUPAC nucleotide ambiguity codes

The single-letter codes for incompletely specified bases were defined
by the Nomenclature Committee of the International Union of Biochemistry
(NC-IUB), Cornish-Bowden (1985), *Eur. J. Biochem.* 150: 1–5
(*Nucl. Acids Res.* 13: 3021–3030):

| Code | Bases | Complement |
|------|-------|------------|
| A | A | T |
| C | C | G |
| G | G | C |
| T | T | A |
| U | U (RNA) | A |
| R | A or G (puRine) | Y |
| Y | C or T (pYrimidine) | R |
| S | G or C (Strong) | S |
| W | A or T (Weak) | W |
| K | G or T (Keto) | M |
| M | A or C (aMino) | K |
| B | C, G or T (not A) | V |
| D | A, G or T (not C) | H |
| H | A, C or T (not G) | D |
| V | A, C or G (not T) | B |
| N | A, C, G or T (aNy) | N |

This table is the complement table used by `reverse_complement()` —
so `revcomp` is correct for ambiguous bases, not just A/C/G/T.
`N` (any base) is the most common ambiguity code in real data and is
the one explicitly exercised by the example file and the tests.

Sources:

- Cornish-Bowden (1985), "Nomenclature for incompletely specified
  bases in nucleic acid sequences", *Nucl. Acids Res.* 13: 3021–3030 —
  https://chemview.epa.gov/chemview/proxy?filename=ncn%2F20211111%2FJ-22-0005%2F1%2F09022526804cf475_REF_+Nomeclature.pdf
- UCSC Genome Browser, "Genome Browser IUPAC Codes",
  https://genome.soe.ucsc.edu/goldenPath/help/iupac.html
- Los Alamos HIV Database, "Codes and Symbols in Sequence Alignments",
  https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html

## 3. GC content

GC% is defined here as `100 * (G + C) / sequence_length`. Ambiguous
bases count in the denominator but not in the numerator — i.e. GC% is
reported over raw sequence length, which is the convention used by
most QC summaries of unmasked sequence. (Some tools exclude `N` from
the denominator; that choice is equally defensible, but it must be
stated. This tool states its choice in the `gc_content()` docstring
and in the README.)

Source for the convention discussion: the NCBI FASTA specification
above (only IUPAC symbols are allowed in sequence data) combined with
standard QC practice; the formula itself is arithmetic, not a cited
claim.

## 4. N50 and L50

N50 is the length of the shortest contig such that contigs of that
length or longer cover at least 50% of the total assembly length; L50
is the number of contigs needed to reach that coverage. Formally:
sort lengths descending, accumulate, and take the length at which the
running sum first reaches half of the total; L50 is the count of
records up to and including that one. N50 is a weighted-median-like
statistic that gives more weight to long sequences than the plain
mean or median.

Worked example (also a unit test): lengths `[10, 20, 30, 40]`, total
100, half 50. Sorted descending `[40, 30, 20, 10]`, cumulative sums
40, 70, 90, 100 — the sum first reaches ≥ 50 at the second record, so
N50 = 30 and L50 = 2.

Note: N50/L50 were designed for genome *assemblies*; applied here to
any FASTA file they are simply length-distribution summaries. The
closely related NG50 uses a known genome size instead of the summed
length — not applicable here, since no reference size is given.

Sources:

- Wikipedia, "N50, L50, and related statistics",
  https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
- UC Davis Bioinformatics Core, "Assessing Assembly Completeness",
  https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/QAQC
- Genomics Aotearoa, "Assembly QC — Long-read Assembly",
  https://genomicsaotearoa.github.io/long-read-assembly-one-day/pages/4_assembly_qc/

## 5. Reverse complement

The reverse complement of a DNA sequence is the sequence read on the
opposite strand: reverse the string, then replace each base by its
complement. It is the correct way to represent the same double-stranded
molecule from the other strand's 5'→3' direction. The implementation
uses the full IUPAC complement table from §2, so ambiguous bases are
handled (e.g. `revcomp("ARYN") = "NRYT"` — verified by hand in the
test suite).

## 6. Prior art (what this tool is not)

- **seqtk** (Heng Li) — the reference C toolkit for FASTA/Q manipulation
  (`seq`, `subseq`, `sample`, `comp`); fast and battle-tested, but a
  compiled binary, not a readable teaching implementation.
  https://github.com/lh3/seqtk
- **Biopython** (`Bio.SeqIO`) — the standard Python library for
  sequence I/O; far broader than this tool (dozens of formats, Entrez,
  alignment). This project intentionally avoids the dependency to stay
  zero-install and readable.
  https://biopython.org/
- **seqkit** (Shen et al. 2016, *PLoS ONE* 11(10): e0163962) — a
  cross-platform Go toolkit with subcommands (`stats`, `grep`, `seq`,
  `subseq`) whose command layout inspired this CLI's subcommand names.
  https://bioinf.shenwei.me/seqkit/

This project re-implements a small, well-defined subset (stats,
filter, extract, GC, reverse complement, conversion) in pure Python
with hand-verified arithmetic, as a learning and utility tool — not as
a replacement for seqtk/seqkit on large files.

## 7. Example data provenance

`examples/sample.fasta` contains five records that are **entirely
synthetic**, generated for this project by repeating short motifs
(`ACGT`, `GGGCCC`, `AAAATTTT`, `ACGTN`, `GATTACA`). No sequence is
claimed to come from any real organism, gene, or accession — the
headers say "synthetic" explicitly. Expected values (lengths 20/30/40/
15/60 bp, GC% 50/100/0/40/28.33, N50 = 40, L50 = 2) were computed by
hand from the motifs and then confirmed by running the tool; the
confirmation is in `examples/example_output.txt`.

## References

1. NCBI. "FASTA Format for Nucleotide Sequences".
   https://www.ncbi.nlm.nih.gov/genbank/fastaformat/
2. NCBI Bookshelf. GenBank Submissions Handbook, "FASTA Formatting".
   https://www.ncbi.nlm.nih.gov/books/NBK566986/
3. Wikipedia. "FASTA format". https://en.wikipedia.org/wiki/FASTA_format
4. Cornish-Bowden, A. (1985). "Nomenclature for incompletely specified
   bases in nucleic acid sequences." *Nucleic Acids Research* 13:
   3021–3030.
   https://chemview.epa.gov/chemview/proxy?filename=ncn%2F20211111%2FJ-22-0005%2F1%2F09022526804cf475_REF_+Nomeclature.pdf
5. UCSC Genome Browser. "Genome Browser IUPAC Codes".
   https://genome.soe.ucsc.edu/goldenPath/help/iupac.html
6. Los Alamos HIV Database. "Codes and Symbols in Sequence Alignments".
   https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html
7. Wikipedia. "N50, L50, and related statistics".
   https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
8. UC Davis Bioinformatics Core. "Assessing Assembly Completeness".
   https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/QAQC
9. Genomics Aotearoa. "Assembly QC — Long-read Assembly".
   https://genomicsaotearoa.github.io/long-read-assembly-one-day/pages/4_assembly_qc/
10. Li, H. seqtk. https://github.com/lh3/seqtk
11. Cock, P. J. A., et al. (2009). "Biopython: freely available Python
    tools for computational molecular biology and bioinformatics."
    *Bioinformatics* 25(11): 1422–1423. https://biopython.org/
12. Shen, W., Le, S., Li, Y., & Hu, F. (2016). "SeqKit: A
    Cross-Platform and Ultrafast Toolkit for FASTA/Q File
    Manipulation." *PLoS ONE* 11(10): e0163962.
    https://bioinf.shenwei.me/seqkit/
