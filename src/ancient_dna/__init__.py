"""
Ancient DNA Analyzer — AI-Powered Paleogenomics Toolkit
========================================================

A comprehensive Python package for analyzing ancient DNA sequences
with genome engineering tools, premium visualizations, and AI-powered insights.

Modules:
    core        — Validation, statistics, FASTA parsing, transcription
    genomics    — Codon analysis, ORF finding, mutation detection,
                  restriction enzymes, protein translation, GC analysis
    visualization — Premium interactive Plotly charts
    ai          — LLM-powered genomic insights via Ollama
    reporting   — Automated report generation
"""

__version__ = "2.0.0"
__author__ = "rudra0812"

from src.ancient_dna.core.validator import validate_dna_sequence
from src.ancient_dna.core.stats import calculate_dna_stats
from src.ancient_dna.core.fasta_parser import read_fasta_file
from src.ancient_dna.core.transcription import (
    transcribe_dna_to_rna,
    reverse_complement,
    complement_strand,
)
