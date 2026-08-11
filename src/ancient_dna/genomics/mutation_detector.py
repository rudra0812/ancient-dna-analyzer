"""
Mutation & SNP Detector
========================

Compares two DNA sequences position-by-position to identify mutations:
single nucleotide polymorphisms (SNPs), and classifies them as transitions
(purine↔purine or pyrimidine↔pyrimidine) or transversions (purine↔pyrimidine).

Critical for paleogenomics — ancient DNA damage patterns typically show
characteristic C→T and G→A transitions caused by cytosine deamination.
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Classification
PURINES = {'A', 'G'}
PYRIMIDINES = {'C', 'T'}


def detect_mutations(seq1: str, seq2: str) -> List[Dict]:
    """
    Detect mutations between two aligned DNA sequences.

    Both sequences must be the same length (pre-aligned).

    Args:
        seq1: Reference DNA sequence.
        seq2: Query DNA sequence (same length as seq1).

    Returns:
        List of mutation dicts with keys:
            - position: 0-indexed position
            - ref_base: Base in seq1
            - alt_base: Base in seq2
            - type: 'transition' or 'transversion'
    """
    s1 = seq1.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    s2 = seq2.upper().replace('\n', '').replace('\r', '').replace(' ', '')

    if len(s1) != len(s2):
        raise ValueError(
            f"Sequences must be the same length for comparison. "
            f"Got {len(s1)} and {len(s2)}."
        )

    mutations: List[Dict] = []
    for i, (b1, b2) in enumerate(zip(s1, s2)):
        if b1 != b2 and b1 in 'ATCG' and b2 in 'ATCG':
            mut_type = _classify_single(b1, b2)
            mutations.append({
                'position': i,
                'ref_base': b1,
                'alt_base': b2,
                'type': mut_type,
            })

    logger.info("Detected %d mutations across %d positions", len(mutations), len(s1))
    return mutations


def _classify_single(base1: str, base2: str) -> str:
    """Classify a single nucleotide change as transition or transversion."""
    if (base1 in PURINES and base2 in PURINES) or \
       (base1 in PYRIMIDINES and base2 in PYRIMIDINES):
        return 'transition'
    return 'transversion'


def mutation_rate(seq1: str, seq2: str) -> Dict[str, float]:
    """
    Calculate mutation rates between two aligned sequences.

    Returns:
        Dict with keys:
            - total_mutations: Count of mismatches
            - mutation_rate: Fraction of positions that differ
            - transitions: Count of transitions
            - transversions: Count of transversions
            - ti_tv_ratio: Transition/transversion ratio
    """
    mutations = detect_mutations(seq1, seq2)
    length = len(seq1.upper().replace('\n', '').replace('\r', '').replace(' ', ''))
    total = len(mutations)

    transitions = sum(1 for m in mutations if m['type'] == 'transition')
    transversions = sum(1 for m in mutations if m['type'] == 'transversion')

    return {
        'total_mutations': total,
        'mutation_rate': round(total / length, 6) if length > 0 else 0.0,
        'transitions': transitions,
        'transversions': transversions,
        'ti_tv_ratio': round(transitions / transversions, 3) if transversions > 0 else float('inf'),
    }


def classify_mutations(seq1: str, seq2: str) -> Dict[str, int]:
    """
    Classify all mutations into specific substitution types.

    Returns:
        Dict mapping substitution pattern (e.g., 'A→G') to count.
        Also includes 'total_transitions' and 'total_transversions'.
    """
    mutations = detect_mutations(seq1, seq2)
    substitutions: Dict[str, int] = {}

    for m in mutations:
        key = f"{m['ref_base']}→{m['alt_base']}"
        substitutions[key] = substitutions.get(key, 0) + 1

    transitions = sum(1 for m in mutations if m['type'] == 'transition')
    transversions = sum(1 for m in mutations if m['type'] == 'transversion')
    substitutions['total_transitions'] = transitions
    substitutions['total_transversions'] = transversions

    return substitutions
