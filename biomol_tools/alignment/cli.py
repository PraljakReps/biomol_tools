"""Command-line interface for alignment tools."""

import argparse
import sys

from .pairwise import global_align, local_align, align_nucleotides


def read_sequence(input_str: str) -> str:
    """Read sequence from file path or return as-is if it's a sequence."""
    if input_str.endswith(('.fasta', '.fa', '.fna', '.faa', '.txt')):
        with open(input_str) as f:
            lines = f.readlines()
            # Skip header lines starting with >
            seq_lines = [l.strip() for l in lines if not l.startswith('>')]
            return ''.join(seq_lines)
    return input_str


def main():
    parser = argparse.ArgumentParser(
        prog='biomol-align',
        description='Pairwise sequence alignment for proteins and nucleotides'
    )
    
    parser.add_argument('seq1', help='First sequence or path to FASTA file')
    parser.add_argument('seq2', help='Second sequence or path to FASTA file')
    
    parser.add_argument(
        '-m', '--mode',
        choices=['global', 'local'],
        default='global',
        help='Alignment mode (default: global)'
    )
    
    parser.add_argument(
        '-t', '--type',
        choices=['protein', 'dna', 'rna'],
        default='protein',
        help='Sequence type (default: protein)'
    )
    
    parser.add_argument(
        '--matrix',
        default='BLOSUM62',
        help='Substitution matrix for proteins (default: BLOSUM62)'
    )
    
    parser.add_argument(
        '--gap-open',
        type=float,
        default=-10.0,
        help='Gap opening penalty (default: -10.0)'
    )
    
    parser.add_argument(
        '--gap-extend',
        type=float,
        default=-0.5,
        help='Gap extension penalty (default: -0.5)'
    )
    
    parser.add_argument(
        '--match',
        type=float,
        default=2.0,
        help='Match score for nucleotides (default: 2.0)'
    )
    
    parser.add_argument(
        '--mismatch',
        type=float,
        default=-1.0,
        help='Mismatch score for nucleotides (default: -1.0)'
    )
    
    parser.add_argument(
        '-s', '--score-only',
        action='store_true',
        help='Only output the alignment score'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output file path (default: print to stdout)'
    )
    
    args = parser.parse_args()
    
    # Read sequences
    seq1 = read_sequence(args.seq1)
    seq2 = read_sequence(args.seq2)
    
    # Perform alignment
    if args.type in ('dna', 'rna'):
        result = align_nucleotides(
            seq1, seq2,
            mode=args.mode,
            match=args.match,
            mismatch=args.mismatch,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
        )
    else:
        if args.mode == 'global':
            result = global_align(
                seq1, seq2,
                matrix=args.matrix,
                gap_open=args.gap_open,
                gap_extend=args.gap_extend,
            )
        else:
            result = local_align(
                seq1, seq2,
                matrix=args.matrix,
                gap_open=args.gap_open,
                gap_extend=args.gap_extend,
            )
    
    # Format output
    if args.score_only:
        output = str(result.score)
    else:
        output = f"{result}\n\nScore: {result.score}"
    
    # Write to file or stdout
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output + '\n')
        print(f"Written to {args.output}")
    else:
        print(output)


if __name__ == '__main__':
    main()
