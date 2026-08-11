"""
DNA Transcription & Strand Operations
======================================

Provides DNA → RNA transcription, reverse complement generation,
and complement strand calculation.
"""

import logging

from src.ancient_dna.core.validator import validate_dna_sequence

logger = logging.getLogger(__name__)

# Complement mapping
_COMPLEMENT_MAP = str.maketrans('ATCG', 'TAGC')


def transcribe_dna_to_rna(dna_sequence: str) -> str:
    """
    Transcribe a DNA sequence into an RNA sequence by replacing T with U.

    Args:
        dna_sequence: A valid DNA sequence string.

    Returns:
        The corresponding RNA sequence.

    Raises:
        ValueError: If the sequence is invalid.
    """
    validate_dna_sequence(dna_sequence)
    rna = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '').replace('T', 'U')
    logger.debug("Transcribed %d bp DNA → RNA", len(rna))
    return rna


def complement_strand(dna_sequence: str) -> str:
    """
    Return the complementary strand of a DNA sequence.

    A↔T, C↔G

    Args:
        dna_sequence: A valid DNA sequence string.

    Returns:
        The complement strand (5'→3' of the complement).
    """
    validate_dna_sequence(dna_sequence)
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    return seq.translate(_COMPLEMENT_MAP)


def reverse_complement(dna_sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    This is the sequence of the complementary strand read in the 3'→5' direction
    (i.e., the antiparallel partner strand read 5'→3').

    Args:
        dna_sequence: A valid DNA sequence string.

    Returns:
        The reverse complement string.
    """
    return complement_strand(dna_sequence)[::-1]
