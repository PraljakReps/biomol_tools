"""
Pairwise sequence alignment utilities for proteins and nucleotides.
"""

from dataclasses import dataclass
from typing import Literal

from Bio import Align
from Bio.Align import substitution_matrices


@dataclass
class AlignmentResult:
    """Container for pairwise alignment results."""
    seq1_aligned: str
    seq2_aligned: str
    score: float
    alignment_str: str
    
    def __repr__(self) -> str:
        return self.alignment_str


def get_aligner(
    mode: Literal["global", "local"] = "global",
    matrix: str = "BLOSUM62",
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> Align.PairwiseAligner:
    """
    Create a configured PairwiseAligner.
    
    Args:
        mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)
        matrix: Substitution matrix name (BLOSUM62, BLOSUM45, PAM250, etc.)
                For nucleotides, use 'NUC.4.4' or set to None for simple match/mismatch
        gap_open: Gap opening penalty (should be negative)
        gap_extend: Gap extension penalty (should be negative)
    
    Returns:
        Configured PairwiseAligner instance
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    
    if matrix:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    return aligner


def align(
    seq1: str,
    seq2: str,
    mode: Literal["global", "local"] = "global",
    matrix: str = "BLOSUM62",
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> AlignmentResult:
    """
    Perform pairwise sequence alignment.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)
        matrix: Substitution matrix name
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
    
    Returns:
        AlignmentResult with aligned sequences and score
    
    Example:
        >>> result = align("MKTLLILAVVAAALA", "MKTLLIFAVVAALA")
        >>> print(result)
        >>> print(f"Score: {result.score}")
    """
    aligner = get_aligner(mode, matrix, gap_open, gap_extend)
    alignments = aligner.align(seq1, seq2)
    
    if not alignments:
        raise ValueError("No alignment found")
    
    best = alignments[0]
    aligned_seqs = best.format().split("\n")
    
    # Parse aligned sequences from the alignment format
    seq1_aligned = str(best).split("\n")[0]
    seq2_aligned = str(best).split("\n")[2]
    
    return AlignmentResult(
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        score=best.score,
        alignment_str=str(best),
    )


def global_align(
    seq1: str,
    seq2: str,
    matrix: str = "BLOSUM62",
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> AlignmentResult:
    """
    Needleman-Wunsch global alignment.
    
    Aligns entire length of both sequences end-to-end.
    Best for sequences of similar length that are expected to be 
    homologous across their full length.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        matrix: Substitution matrix name
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
    
    Returns:
        AlignmentResult
    """
    return align(seq1, seq2, mode="global", matrix=matrix, 
                 gap_open=gap_open, gap_extend=gap_extend)


def local_align(
    seq1: str,
    seq2: str,
    matrix: str = "BLOSUM62",
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> AlignmentResult:
    """
    Smith-Waterman local alignment.
    
    Finds the best matching subsequence between the two sequences.
    Best for finding conserved domains or motifs, or when sequences
    have different lengths or only partial homology.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        matrix: Substitution matrix name
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
    
    Returns:
        AlignmentResult
    """
    return align(seq1, seq2, mode="local", matrix=matrix,
                 gap_open=gap_open, gap_extend=gap_extend)


def align_nucleotides(
    seq1: str,
    seq2: str,
    mode: Literal["global", "local"] = "global",
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -5.0,
    gap_extend: float = -0.5,
) -> AlignmentResult:
    """
    Align DNA/RNA sequences with simple match/mismatch scoring.
    
    Args:
        seq1: First nucleotide sequence
        seq2: Second nucleotide sequence
        mode: 'global' or 'local'
        match: Score for matching bases
        mismatch: Score for mismatching bases
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
    
    Returns:
        AlignmentResult
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    alignments = aligner.align(seq1.upper(), seq2.upper())
    
    if not alignments:
        raise ValueError("No alignment found")
    
    best = alignments[0]
    seq1_aligned = str(best).split("\n")[0]
    seq2_aligned = str(best).split("\n")[2]
    
    return AlignmentResult(
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        score=best.score,
        alignment_str=str(best),
    )


# Convenience aliases
nw_align = global_align  # Needleman-Wunsch
sw_align = local_align   # Smith-Waterman


if __name__ == "__main__":
    # Example usage
    protein1 = "MKTLLILAVVAAALA"
    protein2 = "MKTLLIFAVVAALA"
    
    print("=== Global Alignment (Protein) ===")
    result = global_align(protein1, protein2)
    print(result)
    print(f"Score: {result.score}\n")
    
    print("=== Local Alignment (Protein) ===")
    result = local_align(protein1, protein2)
    print(result)
    print(f"Score: {result.score}\n")
    
    dna1 = "ATGCGATCGATCGATCG"
    dna2 = "ATGCAATCGTTCGATCG"
    
    print("=== Global Alignment (DNA) ===")
    result = align_nucleotides(dna1, dna2, mode="global")
    print(result)
    print(f"Score: {result.score}")
