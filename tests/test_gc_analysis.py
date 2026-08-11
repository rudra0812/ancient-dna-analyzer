"""Tests for GC analysis module."""

import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ancient_dna.genomics.gc_analysis import (
    gc_skew,
    sliding_window_gc,
    cumulative_gc_skew,
)


class TestSlidingWindowGC:
    """Tests for sliding_window_gc."""

    def test_returns_list(self):
        result = sliding_window_gc("ATCGATCGATCGATCG" * 10, window_size=20)
        assert isinstance(result, list)

    def test_all_gc_sequence(self):
        result = sliding_window_gc("GCGCGCGCGCGCGCGCGCGC", window_size=10)
        assert all(d['gc_content'] == 100.0 for d in result)

    def test_no_gc_sequence(self):
        result = sliding_window_gc("ATATATATATATATATATATAT", window_size=10)
        assert all(d['gc_content'] == 0.0 for d in result)

    def test_short_sequence(self):
        result = sliding_window_gc("ATCG", window_size=100)
        assert len(result) == 1
        assert result[0]['gc_content'] == 50.0

    def test_has_position_key(self):
        result = sliding_window_gc("ATCGATCGATCG" * 10, window_size=20)
        assert 'position' in result[0]
        assert 'gc_content' in result[0]


class TestGCSkew:
    """Tests for gc_skew."""

    def test_returns_list(self):
        result = gc_skew("ATCGATCG" * 20, window_size=20)
        assert isinstance(result, list)

    def test_all_g_positive_skew(self):
        result = gc_skew("GGGGGGGGGGGGGGGGGGGG", window_size=10)
        assert all(d['skew_value'] > 0 for d in result)

    def test_all_c_negative_skew(self):
        result = gc_skew("CCCCCCCCCCCCCCCCCCCC", window_size=10)
        assert all(d['skew_value'] < 0 for d in result)

    def test_balanced_zero_skew(self):
        result = gc_skew("GCGCGCGCGCGCGCGCGCGC", window_size=10)
        assert all(abs(d['skew_value']) < 0.01 for d in result)


class TestCumulativeGCSkew:
    """Tests for cumulative_gc_skew."""

    def test_returns_list(self):
        result = cumulative_gc_skew("ATCGATCG")
        assert isinstance(result, list)
        assert len(result) == 8

    def test_all_g_increasing(self):
        result = cumulative_gc_skew("GGGG")
        values = [d['cumulative_skew'] for d in result]
        assert values == [1, 2, 3, 4]

    def test_all_c_decreasing(self):
        result = cumulative_gc_skew("CCCC")
        values = [d['cumulative_skew'] for d in result]
        assert values == [-1, -2, -3, -4]

    def test_a_and_t_no_change(self):
        result = cumulative_gc_skew("ATAT")
        values = [d['cumulative_skew'] for d in result]
        assert values == [0, 0, 0, 0]
