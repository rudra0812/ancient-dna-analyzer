"""
Ancient DNA Analyzer — Command Line Interface
===============================================

Usage:
    python cli.py analyze <fasta_file>        Analyze sequences from a FASTA file
    python cli.py analyze --sequence <seq>     Analyze a raw DNA sequence
    python cli.py compare <file1> <file2>      Compare two sequences
    python cli.py info                         Show tool information
"""

import argparse
import sys
import os

# Fix Windows console encoding for Unicode output
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.ancient_dna.core.validator import validate_dna_sequence
from src.ancient_dna.core.stats import calculate_dna_stats
from src.ancient_dna.core.fasta_parser import read_fasta_file
from src.ancient_dna.core.transcription import transcribe_dna_to_rna, reverse_complement
from src.ancient_dna.genomics.orf_finder import find_orfs, orf_summary
from src.ancient_dna.genomics.codon_analysis import codon_usage_table, codon_bias_score
from src.ancient_dna.genomics.protein_translation import six_frame_translation, find_proteins
from src.ancient_dna.genomics.gc_analysis import sliding_window_gc
from src.ancient_dna.genomics.restriction_enzymes import restriction_map
from src.ancient_dna.genomics.mutation_detector import mutation_rate, classify_mutations
from src.ancient_dna.reporting.report_generator import generate_dna_report


BANNER = """
    +===============================================================+
    |       ANCIENT DNA ANALYZER  v2.0                              |
    |       AI-Powered Paleogenomics Toolkit                        |
    +===============================================================+
"""


def cmd_analyze(args):
    """Run full analysis on sequence(s)."""
    sequences = {}

    if args.sequence:
        sequences['CLI_Input'] = args.sequence.upper()
    elif args.fasta:
        sequences = read_fasta_file(args.fasta)
        if not sequences:
            print("❌ Failed to read FASTA file.")
            sys.exit(1)
    else:
        print("❌ Provide either --fasta or --sequence.")
        sys.exit(1)

    for name, seq in sequences.items():
        print(f"\n{'='*60}")
        print(f"Analyzing: {name} ({len(seq)} bp)")
        print(f"{'='*60}")

        # Basic stats
        stats = calculate_dna_stats(seq)
        if not stats:
            continue

        print(f"  GC Content:    {stats['gc_content']}%")
        print(f"  AT/GC Ratio:   {stats['at_gc_ratio']}")
        print(f"  Pu/Py Ratio:   {stats['purine_pyrimidine_ratio']}")

        # ORFs
        orf_info = orf_summary(seq)
        print(f"  ORFs Found:    {orf_info['total_orfs']} (longest: {orf_info['longest_orf_length']} bp)")

        # Proteins
        proteins = find_proteins(seq, min_length=10)
        print(f"  Proteins:      {len(proteins)} potential protein(s)")

        # Restriction sites
        rmap = restriction_map(seq)
        print(f"  Restriction:   {rmap['total_cut_sites']} sites ({len(rmap['enzymes_found'])} enzymes)")

        # Codon bias
        bias = codon_bias_score(seq)
        print(f"  Codon Bias:    {bias} (20=extreme, 61=none)")

        # Full report
        if args.report:
            generate_dna_report(name, seq, save_report=True, output_dir=args.output or ".")

        print()


def cmd_compare(args):
    """Compare two sequences for mutations."""
    seqs1 = read_fasta_file(args.file1)
    seqs2 = read_fasta_file(args.file2)

    if not seqs1 or not seqs2:
        print("❌ Failed to read one or both FASTA files.")
        sys.exit(1)

    name1, seq1 = list(seqs1.items())[0]
    name2, seq2 = list(seqs2.items())[0]

    min_len = min(len(seq1), len(seq2))
    print(f"\nComparing: {name1} vs {name2}")
    print(f"Alignment length: {min_len} bp")

    rates = mutation_rate(seq1[:min_len], seq2[:min_len])
    print(f"  Total mutations: {rates['total_mutations']}")
    print(f"  Mutation rate:   {rates['mutation_rate']:.4%}")
    print(f"  Transitions:     {rates['transitions']}")
    print(f"  Transversions:   {rates['transversions']}")
    print(f"  Ti/Tv ratio:     {rates['ti_tv_ratio']}")

    subs = classify_mutations(seq1[:min_len], seq2[:min_len])
    print("\n  Substitution types:")
    for key, val in sorted(subs.items()):
        if not key.startswith('total_'):
            print(f"    {key}: {val}")


def cmd_info(args):
    """Show tool information."""
    print(BANNER)
    print("  Modules:")
    print("    Core        — Validation, statistics, FASTA parsing, transcription")
    print("    Genomics    — Codon analysis, ORF finder, mutation detector,")
    print("                  restriction enzymes, protein translation, GC analysis")
    print("    Viz         — 8 premium interactive Plotly charts")
    print("    AI          — LLM-powered insights via Ollama")
    print("    Reporting   — Automated text report generation")
    print()
    print("  GitHub: https://github.com/rudra0812/ancient-dna-analyzer")
    print()


def main():
    parser = argparse.ArgumentParser(
        prog='ancient-dna',
        description='🧬 Ancient DNA Analyzer — AI-Powered Paleogenomics Toolkit',
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # analyze
    p_analyze = subparsers.add_parser('analyze', help='Analyze DNA sequence(s)')
    p_analyze.add_argument('--fasta', '-f', type=str, help='Path to FASTA file')
    p_analyze.add_argument('--sequence', '-s', type=str, help='Raw DNA sequence string')
    p_analyze.add_argument('--report', '-r', action='store_true', help='Generate text report')
    p_analyze.add_argument('--output', '-o', type=str, help='Output directory')
    p_analyze.set_defaults(func=cmd_analyze)

    # compare
    p_compare = subparsers.add_parser('compare', help='Compare two sequences')
    p_compare.add_argument('file1', type=str, help='First FASTA file')
    p_compare.add_argument('file2', type=str, help='Second FASTA file')
    p_compare.set_defaults(func=cmd_compare)

    # info
    p_info = subparsers.add_parser('info', help='Show tool information')
    p_info.set_defaults(func=cmd_info)

    args = parser.parse_args()

    print(BANNER)

    if not args.command:
        parser.print_help()
        return

    args.func(args)


if __name__ == '__main__':
    main()
