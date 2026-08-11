"""
Open Reading Frame (ORF) Finder
================================

Scans DNA sequences across all 6 reading frames (3 forward + 3 reverse complement)
to identify potential protein-coding regions. An ORF is defined as a region
bounded by a start codon (ATG) and a stop codon (TAA, TAG, TGA).

Essential for gene prediction in fragmented ancient DNA where complete gene
annotations may be unavailable.
"""

import logging
from typing import Dict, List, Optional

from src.ancient_dna.core.transcription import reverse_complement

logger = logging.getLogger(__name__)

START_CODON = 'ATG'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def find_orfs(
    dna_sequence: str,
    min_length: int = 30,
    both_strands: bool = True,
) -> List[Dict]:
    """
    Find all Open Reading Frames in a DNA sequence.

    Args:
        dna_sequence: DNA sequence string.
        min_length: Minimum ORF length in nucleotides (default 30 = 10 amino acids).
        both_strands: If True, search both forward and reverse complement strands.

    Returns:
        List of dicts with keys:
            - start: Start position (0-indexed, in the original sequence)
            - end: End position (exclusive)
            - length: ORF length in nucleotides
            - frame: Reading frame (1, 2, 3 for forward; -1, -2, -3 for reverse)
            - strand: '+' or '-'
            - sequence: The ORF nucleotide sequence
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    orfs: List[Dict] = []

    # Search forward strand (frames +1, +2, +3)
    for frame in range(3):
        _scan_frame(seq, frame, frame + 1, '+', min_length, orfs)

    # Search reverse complement (frames -1, -2, -3)
    if both_strands:
        rev_seq = reverse_complement(seq)
        for frame in range(3):
            _scan_frame(rev_seq, frame, -(frame + 1), '-', min_length, orfs)

    # Sort by length descending
    orfs.sort(key=lambda o: o['length'], reverse=True)
    logger.info("Found %d ORFs (min_length=%d bp)", len(orfs), min_length)
    return orfs


def _scan_frame(
    seq: str, frame_offset: int, frame_label: int,
    strand: str, min_length: int, orfs: List[Dict]
) -> None:
    """Scan a single reading frame for ORFs."""
    i = frame_offset
    while i < len(seq) - 2:
        codon = seq[i:i + 3]
        if codon == START_CODON:
            # Found a start codon — scan for stop
            for j in range(i + 3, len(seq) - 2, 3):
                stop = seq[j:j + 3]
                if stop in STOP_CODONS:
                    orf_len = j + 3 - i
                    if orf_len >= min_length:
                        orfs.append({
                            'start': i,
                            'end': j + 3,
                            'length': orf_len,
                            'frame': frame_label,
                            'strand': strand,
                            'sequence': seq[i:j + 3],
                        })
                    i = j + 3  # Continue after this stop codon
                    break
            else:
                # No stop codon found — move to next codon
                i += 3
                continue
        else:
            i += 3


def get_longest_orf(dna_sequence: str, both_strands: bool = True) -> Optional[Dict]:
    """
    Get the longest ORF in a DNA sequence.

    Args:
        dna_sequence: DNA sequence string.
        both_strands: Search both strands.

    Returns:
        Dict for the longest ORF, or None if no ORFs found.
    """
    orfs = find_orfs(dna_sequence, min_length=3, both_strands=both_strands)
    return orfs[0] if orfs else None


def orf_summary(dna_sequence: str, min_length: int = 30) -> Dict:
    """
    Generate a summary of ORFs found in a sequence.

    Returns:
        Dict with keys: total_orfs, forward_orfs, reverse_orfs,
        longest_orf_length, avg_orf_length, orfs_by_frame
    """
    orfs = find_orfs(dna_sequence, min_length=min_length)
    if not orfs:
        return {
            'total_orfs': 0, 'forward_orfs': 0, 'reverse_orfs': 0,
            'longest_orf_length': 0, 'avg_orf_length': 0, 'orfs_by_frame': {},
        }

    forward = [o for o in orfs if o['strand'] == '+']
    reverse = [o for o in orfs if o['strand'] == '-']
    lengths = [o['length'] for o in orfs]

    by_frame: Dict[int, int] = {}
    for o in orfs:
        by_frame[o['frame']] = by_frame.get(o['frame'], 0) + 1

    return {
        'total_orfs': len(orfs),
        'forward_orfs': len(forward),
        'reverse_orfs': len(reverse),
        'longest_orf_length': max(lengths),
        'avg_orf_length': round(sum(lengths) / len(lengths), 1),
        'orfs_by_frame': by_frame,
    }
