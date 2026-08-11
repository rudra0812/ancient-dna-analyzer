"""Tests for codon analysis module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.codon_analysis import (
    count_codons,
    codon_usage_table,
    codon_bias_score,
)


class TestCountCodons:
    """Tests for count_codons."""

    def test_basic_counting(self):
        counts = count_codons("ATGATGATG")
        assert counts['ATG'] == 3

    def test_frame_0(self):
        counts = count_codons("ATGAAACCC", reading_frame=0)
        assert 'ATG' in counts
        assert 'AAA' in counts
        assert 'CCC' in counts

    def test_frame_1(self):
        counts = count_codons("AATGAAACCC", reading_frame=1)
        assert 'ATG' in counts

    def test_invalid_frame_raises(self):
        with pytest.raises(ValueError):
            count_codons("ATGATG", reading_frame=3)

    def test_empty_sequence(self):
        counts = count_codons("")
        assert counts == {}

    def test_short_sequence(self):
        counts = count_codons("AT")
        assert counts == {}


class TestCodonUsageTable:
    """Tests for codon_usage_table."""

    def test_returns_list(self):
        table = codon_usage_table("ATGATGATG")
        assert isinstance(table, list)
        assert len(table) == 64  # All 64 codons

    def test_has_required_keys(self):
        table = codon_usage_table("ATGATGATG")
        entry = table[0]
        assert 'codon' in entry
        assert 'amino_acid' in entry
        assert 'count' in entry
        assert 'frequency' in entry
        assert 'rscu' in entry


class TestCodonBiasScore:
    """Tests for codon_bias_score."""

    def test_returns_float(self):
        score = codon_bias_score("ATGATGATGATGATG")
        assert isinstance(score, float)

    def test_score_range(self):
        score = codon_bias_score("ATGATGATGATGATGAAACCC")
        assert 1.0 <= score <= 61.0

    def test_empty_sequence(self):
        score = codon_bias_score("")
        assert score == 61.0  # No data → no bias
