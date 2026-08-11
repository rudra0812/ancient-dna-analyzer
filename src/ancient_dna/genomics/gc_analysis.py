"""
GC Content & Skew Analysis
============================

Performs sliding-window GC content analysis and GC skew calculations.

GC skew = (G - C) / (G + C) measures the asymmetry of G and C bases.
In bacterial genomes, the cumulative GC skew plot reveals the origin
and terminus of replication. In ancient DNA, unusual GC patterns can
indicate post-mortem damage or contamination.
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def sliding_window_gc(
    dna_sequence: str,
    window_size: int = 100,
    step_size: int = 10,
) -> List[Dict]:
    """
    Calculate GC content using a sliding window across the sequence.

    Args:
        dna_sequence: DNA sequence string.
        window_size: Size of the sliding window in base pairs.
        step_size: Number of bases to slide the window by.

    Returns:
        List of dicts with keys: position (center of window), gc_content
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    results: List[Dict] = []

    if len(seq) < window_size:
        # Sequence shorter than window — compute single value
        gc = _gc_of(seq)
        return [{'position': len(seq) // 2, 'gc_content': gc}]

    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i + window_size]
        gc = _gc_of(window)
        results.append({
            'position': i + window_size // 2,
            'gc_content': gc,
        })

    logger.debug("Sliding window GC: %d windows (size=%d, step=%d)",
                 len(results), window_size, step_size)
    return results


def gc_skew(
    dna_sequence: str,
    window_size: int = 100,
    step_size: int = 10,
) -> List[Dict]:
    """
    Calculate GC skew (G-C)/(G+C) using a sliding window.

    Args:
        dna_sequence: DNA sequence string.
        window_size: Sliding window size.
        step_size: Step size for the window.

    Returns:
        List of dicts with keys: position, skew_value
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    results: List[Dict] = []

    if len(seq) < window_size:
        skew = _skew_of(seq)
        return [{'position': len(seq) // 2, 'skew_value': skew}]

    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i + window_size]
        skew = _skew_of(window)
        results.append({
            'position': i + window_size // 2,
            'skew_value': skew,
        })

    return results


def cumulative_gc_skew(dna_sequence: str) -> List[Dict]:
    """
    Calculate cumulative GC skew across the sequence (per-base resolution).

    The cumulative GC skew is the running sum of (G - C) at each position.
    The minimum of this curve indicates the origin of replication;
    the maximum indicates the terminus.

    Returns:
        List of dicts with keys: position, cumulative_skew
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    results: List[Dict] = []
    cum_skew = 0.0

    for i, base in enumerate(seq):
        if base == 'G':
            cum_skew += 1
        elif base == 'C':
            cum_skew -= 1
        results.append({
            'position': i,
            'cumulative_skew': cum_skew,
        })

    return results


def _gc_of(seq: str) -> float:
    """Calculate GC percentage for a sequence."""
    gc = sum(1 for b in seq if b in 'GC')
    return round((gc / len(seq)) * 100, 2) if len(seq) > 0 else 0.0


def _skew_of(seq: str) -> float:
    """Calculate GC skew for a sequence."""
    g = seq.count('G')
    c = seq.count('C')
    denominator = g + c
    return round((g - c) / denominator, 4) if denominator > 0 else 0.0
