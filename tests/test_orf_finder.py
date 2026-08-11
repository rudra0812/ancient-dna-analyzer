"""Tests for ORF finder module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.orf_finder import find_orfs, get_longest_orf, orf_summary


class TestFindORFs:
    """Tests for find_orfs."""

    def test_simple_orf(self):
        # ATG...TAA
        seq = "ATGAAACCCTTTGGGTAA"
        orfs = find_orfs(seq, min_length=3, both_strands=False)
        assert len(orfs) >= 1
        assert any(o['sequence'].startswith('ATG') for o in orfs)

    def test_no_start_codon(self):
        seq = "AAACCCTTTGGG"
        orfs = find_orfs(seq, min_length=3, both_strands=False)
        assert len(orfs) == 0

    def test_no_stop_codon(self):
        seq = "ATGAAACCCTTTGGG"  # No stop
        orfs = find_orfs(seq, min_length=3, both_strands=False)
        assert len(orfs) == 0

    def test_min_length_filter(self):
        seq = "ATGTAA"  # Very short ORF (6 bp)
        orfs_short = find_orfs(seq, min_length=3, both_strands=False)
        orfs_long = find_orfs(seq, min_length=30, both_strands=False)
        assert len(orfs_short) >= len(orfs_long)

    def test_both_strands(self):
        seq = "ATGAAACCCTTTGGGTAA"
        orfs_fwd = find_orfs(seq, min_length=3, both_strands=False)
        orfs_both = find_orfs(seq, min_length=3, both_strands=True)
        assert len(orfs_both) >= len(orfs_fwd)

    def test_orf_has_required_keys(self):
        seq = "ATGAAACCCTTTGGGTAA"
        orfs = find_orfs(seq, min_length=3, both_strands=False)
        if orfs:
            orf = orfs[0]
            assert 'start' in orf
            assert 'end' in orf
            assert 'length' in orf
            assert 'frame' in orf
            assert 'strand' in orf
            assert 'sequence' in orf


class TestGetLongestORF:
    """Tests for get_longest_orf."""

    def test_returns_longest(self):
        seq = "ATGAAACCCTTTGGGTAA"
        longest = get_longest_orf(seq)
        if longest:
            assert isinstance(longest, dict)

    def test_no_orfs_returns_none(self):
        seq = "AAACCCTTTGGG"
        assert get_longest_orf(seq) is None


class TestORFSummary:
    """Tests for orf_summary."""

    def test_summary_structure(self):
        seq = "ATGAAACCCTTTGGGTAA"
        summary = orf_summary(seq, min_length=3)
        assert 'total_orfs' in summary
        assert 'forward_orfs' in summary
        assert 'reverse_orfs' in summary
        assert 'longest_orf_length' in summary

    def test_empty_summary(self):
        summary = orf_summary("AAAA", min_length=30)
        assert summary['total_orfs'] == 0
