"""Tests for DNA transcription and strand operations."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.core.transcription import (
    transcribe_dna_to_rna,
    complement_strand,
    reverse_complement,
)


class TestTranscription:
    """Tests for transcribe_dna_to_rna."""

    def test_basic_transcription(self):
        assert transcribe_dna_to_rna("ATCG") == "AUCG"

    def test_no_thymine(self):
        assert transcribe_dna_to_rna("ACGACG") == "ACGACG"

    def test_all_thymine(self):
        assert transcribe_dna_to_rna("TTTT") == "UUUU"

    def test_lowercase_input(self):
        assert transcribe_dna_to_rna("atcg") == "AUCG"

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError):
            transcribe_dna_to_rna("ATCGXYZ")


class TestComplementStrand:
    """Tests for complement_strand."""

    def test_basic_complement(self):
        assert complement_strand("ATCG") == "TAGC"

    def test_complement_all_a(self):
        assert complement_strand("AAAA") == "TTTT"

    def test_complement_palindrome(self):
        assert complement_strand("GAATTC") == "CTTAAG"


class TestReverseComplement:
    """Tests for reverse_complement."""

    def test_basic_reverse_complement(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_palindrome(self):
        # EcoRI site GAATTC — reverse complement is also GAATTC
        assert reverse_complement("GAATTC") == "GAATTC"

    def test_single_base(self):
        assert reverse_complement("A") == "T"

    def test_reverse_complement_of_reverse_complement(self):
        seq = "ATCGATCG"
        assert reverse_complement(reverse_complement(seq)) == seq
