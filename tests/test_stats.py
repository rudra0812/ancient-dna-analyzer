"""Tests for DNA statistics calculator."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.core.stats import calculate_dna_stats


class TestCalculateDNAStats:
    """Tests for calculate_dna_stats."""

    def test_basic_counts(self):
        stats = calculate_dna_stats("AATTCCGG")
        assert stats['counts'] == {'A': 2, 'T': 2, 'C': 2, 'G': 2}

    def test_length(self):
        stats = calculate_dna_stats("ATCGATCG")
        assert stats['length'] == 8

    def test_percentages_equal_distribution(self):
        stats = calculate_dna_stats("ATCG")
        assert stats['percentages'] == {'A': 25.0, 'T': 25.0, 'C': 25.0, 'G': 25.0}

    def test_gc_content_50(self):
        stats = calculate_dna_stats("ATCG")
        assert stats['gc_content'] == 50.0

    def test_gc_content_100(self):
        stats = calculate_dna_stats("CCGG")
        assert stats['gc_content'] == 100.0

    def test_gc_content_0(self):
        stats = calculate_dna_stats("AATT")
        assert stats['gc_content'] == 0.0

    def test_at_gc_ratio(self):
        stats = calculate_dna_stats("AATTCCGG")
        assert stats['at_gc_ratio'] == 1.0

    def test_purine_pyrimidine_ratio(self):
        stats = calculate_dna_stats("AAGG" + "CCTT")
        assert stats['purine_pyrimidine_ratio'] == 1.0

    def test_dinucleotide_frequencies(self):
        stats = calculate_dna_stats("ATCG")
        assert 'dinucleotide_freq' in stats
        assert 'AT' in stats['dinucleotide_freq']

    def test_empty_sequence_returns_none(self):
        result = calculate_dna_stats("")
        assert result is None

    def test_lowercase_input(self):
        stats = calculate_dna_stats("atcg")
        assert stats['length'] == 4
        assert stats['counts']['A'] == 1

    def test_single_base_sequence(self):
        stats = calculate_dna_stats("A")
        assert stats['counts']['A'] == 1
        assert stats['gc_content'] == 0.0

    def test_all_same_base(self):
        stats = calculate_dna_stats("AAAAAAAAAA")
        assert stats['percentages']['A'] == 100.0
        assert stats['gc_content'] == 0.0

    def test_whitespace_handled(self):
        stats = calculate_dna_stats("AT CG\nAT CG")
        assert stats['length'] == 8
