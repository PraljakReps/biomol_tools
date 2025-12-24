# biomol_tools

Simple bioinformatics utilities for protein, RNA, and DNA sequence analysis.

## Installation

```bash
pip install biomol_tools
```

Or install from source:

```bash
git clone https://github.com/nprber/biomol_tools.git
cd biomol_tools
pip install -e .
```

## Quick Start

### Python API

```python
from biomol_tools.alignment import global_align, local_align, align_nucleotides

# Global alignment (Needleman-Wunsch)
result = global_align("MKTLLILAVVAAALA", "MKTLLIFAVVAALA")
print(result)
print(f"Score: {result.score}")

# Local alignment (Smith-Waterman)
result = local_align("MKTLLILAVVAAALA", "MKTLLIFAVVAALA")
print(result)

# DNA/RNA alignment
result = align_nucleotides("ATGCGATCGATCG", "ATGCAATCGTTCG", mode="global")
print(result)
```

### Command Line

```bash
# Protein alignment (global by default)
biomol-align MKTLLILAVVAAALA MKTLLIFAVVAALA

# Local alignment
biomol-align MKTLLILAVVAAALA MKTLLIFAVVAALA -m local

# DNA alignment
biomol-align ATGCGATCGATCG ATGCAATCGTTCG -t dna

# From FASTA files
biomol-align seq1.fasta seq2.fasta

# Just the score
biomol-align MKTLLILAVVAAALA MKTLLIFAVVAALA -s

# Custom parameters
biomol-align SEQ1 SEQ2 --matrix BLOSUM45 --gap-open -12 --gap-extend -1
```

Run `biomol-align --help` for all options.

## Features

### Alignment

| Function | Description |
|----------|-------------|
| `global_align()` | Needleman-Wunsch global alignment |
| `local_align()` | Smith-Waterman local alignment |
| `align_nucleotides()` | DNA/RNA alignment with match/mismatch scoring |

**Aliases:** `nw_align()` for global, `sw_align()` for local.

### Substitution Matrices

For proteins, any Biopython-supported matrix works: `BLOSUM62` (default), `BLOSUM45`, `BLOSUM80`, `PAM250`, etc.

## Requirements

- Python ≥ 3.9
- Biopython ≥ 1.80

## License

MIT
