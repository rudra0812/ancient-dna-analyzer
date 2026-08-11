"""Tests for mutation detector module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.mutation_detector import (
    detect_mutations,
    mutation_rate,
    classify_mutations,
)


class TestDetectMutations:
    """Tests for detect_mutations."""

    def test_identical_sequences(self):
        mutations = detect_mutations("ATCGATCG", "ATCGATCG")
        assert len(mutations) == 0

    def test_single_mutation(self):
        mutations = detect_mutations("ATCGATCG", "ATCAATCG")
        assert len(mutations) == 1
        assert mutations[0]['position'] == 3
        assert mutations[0]['ref_base'] == 'G'
        assert mutations[0]['alt_base'] == 'A'

    def test_multiple_mutations(self):
        mutations = detect_mutations("AAAA", "TTTT")
        assert len(mutations) == 4

    def test_transition_classified(self):
        # A → G is a transition (purine → purine)
        mutations = detect_mutations("ATCG", "GTCG")
        assert mutations[0]['type'] == 'transition'

    def test_transversion_classified(self):
        # A → C is a transversion (purine → pyrimidine)
        mutations = detect_mutations("ATCG", "CTCG")
        assert mutations[0]['type'] == 'transversion'

    def test_different_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            detect_mutations("ATCG", "ATC")


class TestMutationRate:
    """Tests for mutation_rate."""

    def test_zero_rate(self):
        rates = mutation_rate("ATCG", "ATCG")
        assert rates['total_mutations'] == 0
        assert rates['mutation_rate'] == 0.0

    def test_100_percent_rate(self):
        rates = mutation_rate("AAAA", "TTTT")
        assert rates['mutation_rate'] == 1.0

    def test_ti_tv_ratio(self):
        # A→G transition, A→C transversion
        rates = mutation_rate("AA", "GC")
        assert rates['transitions'] == 1
        assert rates['transversions'] == 1
        assert rates['ti_tv_ratio'] == 1.0


class TestClassifyMutations:
    """Tests for classify_mutations."""

    def test_substitution_types(self):
        subs = classify_mutations("AG", "GA")
        assert "A→G" in subs
        assert "G→A" in subs

    def test_total_counts(self):
        subs = classify_mutations("AATT", "GGCC")
        assert subs['total_transitions'] + subs['total_transversions'] == 4
