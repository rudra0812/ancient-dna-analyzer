"""
Codon Usage Analysis
====================

Analyzes codon usage patterns in DNA sequences. Codons are three-nucleotide
sequences that encode amino acids. Different organisms show distinct codon
usage biases — analyzing these patterns can reveal evolutionary relationships,
gene expression levels, and help optimize synthetic gene design.

Key concepts:
    - Codon: A triplet of nucleotides (e.g., ATG, GCA) that codes for an amino acid
    - Codon bias: The preferential use of certain codons over synonymous alternatives
    - RSCU (Relative Synonymous Codon Usage): Measures deviation from uniform usage
"""

import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Standard genetic code: codon → amino acid
CODON_TABLE: Dict[str, str] = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stop', 'TAG': 'Stop',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stop', 'TGG': 'Trp',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

# Reverse mapping: amino acid → list of synonymous codons
_AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def count_codons(dna_sequence: str, reading_frame: int = 0) -> Dict[str, int]:
    """
    Count the occurrences of each codon in a DNA sequence.

    Args:
        dna_sequence: DNA sequence string.
        reading_frame: Reading frame offset (0, 1, or 2).

    Returns:
        Dictionary mapping codon → count.
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    if reading_frame not in (0, 1, 2):
        raise ValueError(f"Reading frame must be 0, 1, or 2, got {reading_frame}")

    counts: Dict[str, int] = {}
    for i in range(reading_frame, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and all(b in 'ATCG' for b in codon):
            counts[codon] = counts.get(codon, 0) + 1

    logger.debug("Counted %d unique codons in frame %d", len(counts), reading_frame)
    return counts


def codon_usage_table(dna_sequence: str, reading_frame: int = 0) -> List[Dict]:
    """
    Generate a codon usage table with frequencies and amino acid mappings.

    Args:
        dna_sequence: DNA sequence string.
        reading_frame: Reading frame offset (0, 1, or 2).

    Returns:
        List of dicts with keys: codon, amino_acid, count, frequency, rscu
    """
    counts = count_codons(dna_sequence, reading_frame)
    total_codons = sum(counts.values()) or 1

    # Calculate RSCU (Relative Synonymous Codon Usage)
    # RSCU = observed / expected, where expected = total_for_aa / num_synonymous
    aa_totals: Dict[str, int] = {}
    for codon, count in counts.items():
        aa = CODON_TABLE.get(codon, '???')
        aa_totals[aa] = aa_totals.get(aa, 0) + count

    table: List[Dict] = []
    for codon in sorted(CODON_TABLE.keys()):
        aa = CODON_TABLE[codon]
        count = counts.get(codon, 0)
        frequency = round(count / total_codons, 4)
        num_synonymous = len(_AA_TO_CODONS.get(aa, [codon]))
        aa_total = aa_totals.get(aa, 0)
        expected = aa_total / num_synonymous if num_synonymous > 0 else 0
        rscu = round(count / expected, 3) if expected > 0 else 0.0

        table.append({
            'codon': codon,
            'amino_acid': aa,
            'count': count,
            'frequency': frequency,
            'rscu': rscu,
        })

    return table


def codon_bias_score(dna_sequence: str, reading_frame: int = 0) -> float:
    """
    Calculate a simple codon bias score (Effective Number of Codons approximation).

    A score closer to 20 indicates extreme bias (few codons used),
    a score near 61 indicates uniform/no bias.

    Args:
        dna_sequence: DNA sequence string.
        reading_frame: Reading frame offset (0, 1, or 2).

    Returns:
        Codon bias score (float).
    """
    counts = count_codons(dna_sequence, reading_frame)
    if not counts:
        return 61.0  # No data → assume no bias

    # Group codons by amino acid
    aa_groups: Dict[str, Dict[str, int]] = {}
    for codon, count in counts.items():
        aa = CODON_TABLE.get(codon, '???')
        if aa == 'Stop' or aa == '???':
            continue
        aa_groups.setdefault(aa, {})[codon] = count

    # Calculate homozygosity (F) for each amino acid family
    f_values: List[float] = []
    for aa, codons in aa_groups.items():
        n = sum(codons.values())
        if n <= 1 or len(_AA_TO_CODONS.get(aa, [])) <= 1:
            continue
        f = sum((c / n) ** 2 for c in codons.values())
        f_values.append(f)

    if not f_values:
        return 61.0

    avg_f = sum(f_values) / len(f_values)
    enc = 1.0 / avg_f if avg_f > 0 else 61.0
    return round(min(enc, 61.0), 2)
