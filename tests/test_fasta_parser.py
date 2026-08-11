"""Tests for FASTA file parser."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.core.fasta_parser import read_fasta_file


class TestReadFastaFile:
    """Tests for read_fasta_file."""

    def test_read_sample_file(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Seq1\nATCGATCG\n>Seq2\nGGGGCCCC\n")
        result = read_fasta_file(str(fasta))
        assert result is not None
        assert len(result) == 2
        assert "Seq1" in result
        assert result["Seq1"] == "ATCGATCG"

    def test_multiline_sequence(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Seq1\nATCG\nATCG\nATCG\n")
        result = read_fasta_file(str(fasta))
        assert result["Seq1"] == "ATCGATCGATCG"

    def test_single_sequence(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">MySeq\nAAAA\n")
        result = read_fasta_file(str(fasta))
        assert len(result) == 1
        assert result["MySeq"] == "AAAA"

    def test_file_not_found(self):
        result = read_fasta_file("nonexistent_file.fasta")
        assert result is None

    def test_empty_file(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        result = read_fasta_file(str(fasta))
        assert result == {}

    def test_lowercase_converted_to_upper(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Seq1\natcgatcg\n")
        result = read_fasta_file(str(fasta))
        assert result["Seq1"] == "ATCGATCG"

    def test_header_with_description(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Seq1 This is a description\nATCG\n")
        result = read_fasta_file(str(fasta))
        assert "Seq1 This is a description" in result

    def test_blank_lines_ignored(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">Seq1\n\nATCG\n\nGGGG\n\n")
        result = read_fasta_file(str(fasta))
        assert result["Seq1"] == "ATCGGGGG"
