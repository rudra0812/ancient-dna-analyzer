"""
DNA Statistics Calculator
=========================

Computes comprehensive nucleotide statistics including base counts,
percentages, GC content, AT/GC ratio, purine/pyrimidine ratio,
and dinucleotide frequencies.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def calculate_dna_stats(dna_sequence: str) -> Optional[Dict]:
    """
    Calculate comprehensive statistics for a DNA sequence.

    Args:
        dna_sequence: String containing DNA sequence (A, T, C, G).

    Returns:
        Dictionary containing:
            - length: Sequence length in base pairs
            - counts: Dict of base counts {A, T, C, G}
            - percentages: Dict of base percentages
            - gc_content: GC content percentage
            - at_gc_ratio: AT/GC ratio
            - purine_pyrimidine_ratio: Purine(A+G) / Pyrimidine(C+T) ratio
            - dinucleotide_freq: Dict of all 16 dinucleotide frequencies
        Returns None if the sequence is empty.
    """
    dna_sequence = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    length = len(dna_sequence)

    if length == 0:
        logger.error("Cannot calculate stats for empty sequence.")
        return None

    # --- Base counts ---
    counts: Dict[str, int] = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in dna_sequence:
        if base in counts:
            counts[base] += 1
        else:
            logger.warning("Non-standard base '%s' ignored in statistics.", base)

    # --- Percentages ---
    percentages = {base: round((count / length) * 100, 2) for base, count in counts.items()}

    # --- GC content ---
    gc_content = round(((counts['G'] + counts['C']) / length) * 100, 2)

    # --- AT/GC ratio ---
    gc_total = counts['G'] + counts['C']
    at_total = counts['A'] + counts['T']
    at_gc_ratio = round(at_total / gc_total, 3) if gc_total > 0 else float('inf')

    # --- Purine / Pyrimidine ratio ---
    purine = counts['A'] + counts['G']
    pyrimidine = counts['C'] + counts['T']
    purine_pyrimidine_ratio = round(purine / pyrimidine, 3) if pyrimidine > 0 else float('inf')

    # --- Dinucleotide frequencies ---
    dinucleotides: Dict[str, int] = {}
    for i in range(length - 1):
        di = dna_sequence[i:i + 2]
        if len(di) == 2 and all(b in 'ATCG' for b in di):
            dinucleotides[di] = dinucleotides.get(di, 0) + 1

    total_di = sum(dinucleotides.values()) if dinucleotides else 1
    dinucleotide_freq = {
        di: round((count / total_di) * 100, 2) for di, count in sorted(dinucleotides.items())
    }

    logger.info("Stats calculated: %d bp, GC=%.2f%%", length, gc_content)
    return {
        'length': length,
        'counts': counts,
        'percentages': percentages,
        'gc_content': gc_content,
        'at_gc_ratio': at_gc_ratio,
        'purine_pyrimidine_ratio': purine_pyrimidine_ratio,
        'dinucleotide_freq': dinucleotide_freq,
    }
