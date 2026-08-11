"""Tests for protein translation module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.protein_translation import (
    translate,
    six_frame_translation,
    find_proteins,
)


class TestTranslate:
    """Tests for translate."""

    def test_start_codon(self):
        protein = translate("ATG")
        assert protein == "M"  # Methionine

    def test_stop_codon(self):
        protein = translate("TAA")
        assert protein == "*"

    def test_known_sequence(self):
        # ATG GAA TAA = Met Glu Stop
        protein = translate("ATGGAATAA")
        assert protein == "ME*"

    def test_reading_frame_1(self):
        protein = translate("AATG", reading_frame=1)
        assert protein == "M"

    def test_invalid_frame_raises(self):
        with pytest.raises(ValueError):
            translate("ATG", reading_frame=5)

    def test_empty_sequence(self):
        protein = translate("")
        assert protein == ""

    def test_incomplete_codon_ignored(self):
        protein = translate("ATGA")  # Last 'A' is incomplete
        assert protein == "M"


class TestSixFrameTranslation:
    """Tests for six_frame_translation."""

    def test_returns_six_frames(self):
        result = six_frame_translation("ATCGATCGATCG")
        assert len(result) == 6
        assert '+1' in result
        assert '-3' in result

    def test_all_frames_are_strings(self):
        result = six_frame_translation("ATCGATCGATCG")
        for frame, protein in result.items():
            assert isinstance(protein, str)


class TestFindProteins:
    """Tests for find_proteins."""

    def test_finds_protein(self):
        # ATG + 10 codons + TAA
        seq = "ATG" + "AAA" * 12 + "TAA"
        proteins = find_proteins(seq, min_length=5)
        assert len(proteins) >= 1

    def test_min_length_filter(self):
        seq = "ATGAAATAA"  # M K * = 2 aa protein
        short = find_proteins(seq, min_length=2)
        long = find_proteins(seq, min_length=10)
        assert len(short) >= len(long)

    def test_protein_has_keys(self):
        seq = "ATG" + "AAA" * 12 + "TAA"
        proteins = find_proteins(seq, min_length=5)
        if proteins:
            p = proteins[0]
            assert 'frame' in p
            assert 'start_aa' in p
            assert 'length' in p
            assert 'sequence' in p
