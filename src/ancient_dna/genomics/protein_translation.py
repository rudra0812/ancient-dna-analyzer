"""
DNA to Protein Translation
===========================

Translates DNA sequences into amino acid chains using the standard
genetic code. Supports all 6 reading frames and can identify
protein-coding sequences from start to stop codons.
"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Standard genetic code: codon → single-letter amino acid
GENETIC_CODE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Full amino acid names for display
AA_NAMES: Dict[str, str] = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'E': 'Glutamic acid', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    '*': 'Stop',
}


def translate(dna_sequence: str, reading_frame: int = 0) -> str:
    """
    Translate a DNA sequence into a protein (amino acid) sequence.

    Args:
        dna_sequence: DNA sequence string.
        reading_frame: Reading frame offset (0, 1, or 2).

    Returns:
        Amino acid sequence string (single-letter codes).
        Stop codons are represented as '*'.
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    if reading_frame not in (0, 1, 2):
        raise ValueError(f"Reading frame must be 0, 1, or 2, got {reading_frame}")

    protein: List[str] = []
    for i in range(reading_frame, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3:
            aa = GENETIC_CODE.get(codon, 'X')  # X for unknown
            protein.append(aa)

    return ''.join(protein)


def six_frame_translation(dna_sequence: str) -> Dict[str, str]:
    """
    Translate a DNA sequence in all 6 reading frames.

    Frames +1, +2, +3 are on the forward strand.
    Frames -1, -2, -3 are on the reverse complement strand.

    Returns:
        Dict mapping frame label ('+1', '+2', '+3', '-1', '-2', '-3')
        to protein sequence string.
    """
    from src.ancient_dna.core.transcription import reverse_complement

    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    rev_seq = reverse_complement(seq)

    return {
        '+1': translate(seq, 0),
        '+2': translate(seq, 1),
        '+3': translate(seq, 2),
        '-1': translate(rev_seq, 0),
        '-2': translate(rev_seq, 1),
        '-3': translate(rev_seq, 2),
    }


def find_proteins(dna_sequence: str, min_length: int = 10) -> List[Dict]:
    """
    Find potential proteins (M...* sequences) across all 6 reading frames.

    Args:
        dna_sequence: DNA sequence string.
        min_length: Minimum protein length in amino acids.

    Returns:
        List of dicts with keys:
            - frame: Reading frame label
            - start_aa: Start position in amino acid numbering
            - length: Protein length in amino acids
            - sequence: Protein sequence string
    """
    translations = six_frame_translation(dna_sequence)
    proteins: List[Dict] = []

    for frame, protein_seq in translations.items():
        # Find all M...* patterns
        i = 0
        while i < len(protein_seq):
            if protein_seq[i] == 'M':
                # Find next stop codon
                stop = protein_seq.find('*', i)
                if stop == -1:
                    # No stop found — take to end
                    prot = protein_seq[i:]
                    if len(prot) >= min_length:
                        proteins.append({
                            'frame': frame,
                            'start_aa': i,
                            'length': len(prot),
                            'sequence': prot,
                        })
                    break
                else:
                    prot = protein_seq[i:stop]
                    if len(prot) >= min_length:
                        proteins.append({
                            'frame': frame,
                            'start_aa': i,
                            'length': len(prot),
                            'sequence': prot,
                        })
                    i = stop + 1
            else:
                i += 1

    proteins.sort(key=lambda p: p['length'], reverse=True)
    logger.info("Found %d proteins (min_length=%d aa)", len(proteins), min_length)
    return proteins
