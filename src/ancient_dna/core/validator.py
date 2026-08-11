"""
DNA Sequence Validator
======================

Validates DNA sequences against standard and IUPAC nucleotide codes.
Supports both strict mode (A, T, C, G only) and relaxed mode (allows
ambiguous IUPAC bases like N, R, Y, etc.).
"""

import logging
from typing import Set

logger = logging.getLogger(__name__)

# Standard DNA bases
VALID_BASES: Set[str] = {'A', 'T', 'C', 'G'}

# Extended IUPAC ambiguity codes for DNA
IUPAC_BASES: Set[str] = {
    'A', 'T', 'C', 'G',    # Standard
    'R',                     # A or G (purine)
    'Y',                     # C or T (pyrimidine)
    'S',                     # G or C (strong)
    'W',                     # A or T (weak)
    'K',                     # G or T (keto)
    'M',                     # A or C (amino)
    'B',                     # C, G, or T (not A)
    'D',                     # A, G, or T (not C)
    'H',                     # A, C, or T (not G)
    'V',                     # A, C, or G (not T)
    'N',                     # Any base
}


def validate_dna_sequence(dna_sequence: str, strict: bool = True) -> bool:
    """
    Validate if a DNA sequence contains only valid nucleotide characters.

    Args:
        dna_sequence: String containing DNA sequence.
        strict: If True, only allows A, T, C, G.
                If False, also allows IUPAC ambiguity codes (N, R, Y, etc.).

    Returns:
        True if the sequence is valid.

    Raises:
        ValueError: If the sequence is empty or contains invalid characters.
    """
    if not dna_sequence or not dna_sequence.strip():
        raise ValueError("DNA sequence cannot be empty.")

    dna_sequence = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    allowed = VALID_BASES if strict else IUPAC_BASES
    seq_chars = set(dna_sequence)

    if not seq_chars.issubset(allowed):
        invalid_chars = seq_chars - allowed
        mode = "strict" if strict else "IUPAC"
        raise ValueError(
            f"Invalid character(s) in {mode} mode: {', '.join(sorted(invalid_chars))}. "
            f"Allowed: {', '.join(sorted(allowed))}"
        )

    logger.debug("Sequence validated (%d bp, mode=%s)", len(dna_sequence),
                 "strict" if strict else "IUPAC")
    return True
