"""Tests for restriction enzyme module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.restriction_enzymes import (
    find_cut_sites,
    digest_sequence,
    restriction_map,
)


class TestFindCutSites:
    """Tests for find_cut_sites."""

    def test_ecori_site(self):
        # EcoRI recognizes GAATTC
        seq = "AAAGAATTCAAA"
        sites = find_cut_sites(seq, ['EcoRI'])
        assert 'EcoRI' in sites
        assert len(sites['EcoRI']) == 1

    def test_no_sites(self):
        seq = "AAAAAAAAAA"
        sites = find_cut_sites(seq, ['EcoRI'])
        assert len(sites) == 0

    def test_multiple_sites(self):
        seq = "GAATTCAAAGAATTC"
        sites = find_cut_sites(seq, ['EcoRI'])
        assert len(sites['EcoRI']) == 2

    def test_unknown_enzyme_warning(self):
        sites = find_cut_sites("ATCG", ['FakeEnzyme'])
        assert len(sites) == 0

    def test_all_enzymes_scan(self):
        # Should not crash when scanning with all enzymes
        seq = "ATCGATCGATCG" * 10
        sites = find_cut_sites(seq)
        assert isinstance(sites, dict)


class TestDigestSequence:
    """Tests for digest_sequence."""

    def test_no_cut_returns_whole(self):
        fragments = digest_sequence("AAAAAAAAAA", "EcoRI")
        assert len(fragments) == 1
        assert fragments[0]['length'] == 10

    def test_single_cut_two_fragments(self):
        seq = "AAAGAATTCAAA"
        fragments = digest_sequence(seq, "EcoRI")
        assert len(fragments) == 2

    def test_fragment_has_keys(self):
        fragments = digest_sequence("AAAGAATTCAAA", "EcoRI")
        f = fragments[0]
        assert 'fragment_number' in f
        assert 'start' in f
        assert 'end' in f
        assert 'length' in f


class TestRestrictionMap:
    """Tests for restriction_map."""

    def test_map_structure(self):
        rmap = restriction_map("GAATTCAAAGGATCC")
        assert 'sequence_length' in rmap
        assert 'enzymes_found' in rmap
        assert 'total_cut_sites' in rmap
        assert 'sites' in rmap
        assert 'fragment_counts' in rmap

    def test_map_with_known_sites(self):
        # Contains EcoRI (GAATTC) and BamHI (GGATCC)
        rmap = restriction_map("GAATTCAAAGGATCC")
        assert rmap['total_cut_sites'] >= 2
