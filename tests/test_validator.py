"""Tests for DNA sequence validator."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.core.validator import validate_dna_sequence


class TestValidateDNASequence:
    """Tests for validate_dna_sequence."""

    def test_valid_sequence(self):
        assert validate_dna_sequence("ATCGATCG") is True

    def test_valid_lowercase(self):
        assert validate_dna_sequence("atcgatcg") is True

    def test_valid_mixed_case(self):
        assert validate_dna_sequence("AtCgAtCg") is True

    def test_valid_long_sequence(self):
        seq = "ATCG" * 1000
        assert validate_dna_sequence(seq) is True

    def test_valid_single_base(self):
        for base in ['A', 'T', 'C', 'G']:
            assert validate_dna_sequence(base) is True

    def test_invalid_characters(self):
        with pytest.raises(ValueError, match="Invalid character"):
            validate_dna_sequence("ATCGXYZ")

    def test_empty_sequence(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            validate_dna_sequence("")

    def test_whitespace_only(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            validate_dna_sequence("   ")

    def test_sequence_with_spaces(self):
        """Spaces should be stripped and sequence validated."""
        assert validate_dna_sequence("ATC GAT CG") is True

    def test_sequence_with_newlines(self):
        assert validate_dna_sequence("ATCG\nATCG") is True

    def test_iupac_strict_mode_rejects(self):
        with pytest.raises(ValueError):
            validate_dna_sequence("ATCGN", strict=True)

    def test_iupac_relaxed_mode_accepts(self):
        assert validate_dna_sequence("ATCGNRYSWKM", strict=False) is True

    def test_numbers_rejected(self):
        with pytest.raises(ValueError):
            validate_dna_sequence("ATCG123")

    def test_rna_base_rejected_in_strict(self):
        with pytest.raises(ValueError):
            validate_dna_sequence("AUCG", strict=True)
