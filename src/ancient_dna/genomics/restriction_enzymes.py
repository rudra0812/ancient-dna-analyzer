"""
Restriction Enzyme Cut-Site Finder
====================================

Searches DNA sequences for recognition sites of common restriction
endonucleases. Returns cut positions and predicted fragment sizes.

Restriction enzymes are molecular scissors that cut DNA at specific
sequences — essential for molecular cloning, DNA fingerprinting,
and mapping ancient DNA fragments.
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Common restriction enzymes: name → (recognition_site, cut_offset_from_start)
# Cut offset is on the top strand, measured from the 5' end of the recognition site
ENZYME_DATABASE: Dict[str, Dict] = {
    'EcoRI':   {'site': 'GAATTC', 'cut': 1, 'description': 'Most common lab enzyme'},
    'BamHI':   {'site': 'GGATCC', 'cut': 1, 'description': 'Used in cloning vectors'},
    'HindIII': {'site': 'AAGCTT', 'cut': 1, 'description': 'Common 6-cutter'},
    'NotI':    {'site': 'GCGGCCGC', 'cut': 2, 'description': 'Rare 8-cutter'},
    'XhoI':    {'site': 'CTCGAG', 'cut': 1, 'description': 'Common cloning enzyme'},
    'SalI':    {'site': 'GTCGAC', 'cut': 1, 'description': 'Compatible with XhoI'},
    'PstI':    {'site': 'CTGCAG', 'cut': 5, 'description': '3\' overhang cutter'},
    'SmaI':    {'site': 'CCCGGG', 'cut': 3, 'description': 'Blunt end cutter'},
    'KpnI':    {'site': 'GGTACC', 'cut': 5, 'description': '3\' overhang cutter'},
    'SacI':    {'site': 'GAGCTC', 'cut': 5, 'description': '3\' overhang cutter'},
    'NcoI':    {'site': 'CCATGG', 'cut': 1, 'description': 'Used for protein expression'},
    'NdeI':    {'site': 'CATATG', 'cut': 2, 'description': 'Contains ATG start codon'},
    'XbaI':    {'site': 'TCTAGA', 'cut': 1, 'description': 'Common cloning enzyme'},
    'SpeI':    {'site': 'ACTAGT', 'cut': 1, 'description': 'Compatible with XbaI'},
    'ApaI':    {'site': 'GGGCCC', 'cut': 5, 'description': '3\' overhang cutter'},
    'MfeI':    {'site': 'CAATTG', 'cut': 1, 'description': 'Compatible with EcoRI'},
}


def find_cut_sites(
    dna_sequence: str,
    enzymes: List[str] = None,
) -> Dict[str, List[int]]:
    """
    Find restriction enzyme cut sites in a DNA sequence.

    Args:
        dna_sequence: DNA sequence string.
        enzymes: List of enzyme names to search for.
                 If None, searches all enzymes in the database.

    Returns:
        Dict mapping enzyme_name → list of cut positions (0-indexed).
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    enzyme_list = enzymes or list(ENZYME_DATABASE.keys())
    results: Dict[str, List[int]] = {}

    for enzyme_name in enzyme_list:
        if enzyme_name not in ENZYME_DATABASE:
            logger.warning("Unknown enzyme: '%s'", enzyme_name)
            continue

        site = ENZYME_DATABASE[enzyme_name]['site']
        positions: List[int] = []

        # Scan for recognition sites
        start = 0
        while True:
            pos = seq.find(site, start)
            if pos == -1:
                break
            cut_pos = pos + ENZYME_DATABASE[enzyme_name]['cut']
            positions.append(cut_pos)
            start = pos + 1

        if positions:
            results[enzyme_name] = positions

    total = sum(len(v) for v in results.values())
    logger.info("Found %d cut sites across %d enzymes", total, len(results))
    return results


def digest_sequence(dna_sequence: str, enzyme_name: str) -> List[Dict]:
    """
    Simulate restriction enzyme digestion and return fragments.

    Args:
        dna_sequence: DNA sequence string.
        enzyme_name: Name of the restriction enzyme.

    Returns:
        List of dicts with keys: fragment_number, start, end, length, sequence_preview
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    sites = find_cut_sites(seq, [enzyme_name])
    cut_positions = sites.get(enzyme_name, [])

    if not cut_positions:
        return [{
            'fragment_number': 1,
            'start': 0,
            'end': len(seq),
            'length': len(seq),
            'sequence_preview': seq[:50] + ('...' if len(seq) > 50 else ''),
        }]

    # Add sequence boundaries
    boundaries = [0] + sorted(cut_positions) + [len(seq)]
    fragments: List[Dict] = []

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]
        frag_seq = seq[start:end]
        fragments.append({
            'fragment_number': i + 1,
            'start': start,
            'end': end,
            'length': end - start,
            'sequence_preview': frag_seq[:50] + ('...' if len(frag_seq) > 50 else ''),
        })

    fragments.sort(key=lambda f: f['length'], reverse=True)
    return fragments


def restriction_map(dna_sequence: str, enzymes: List[str] = None) -> Dict:
    """
    Generate a complete restriction map for a DNA sequence.

    Returns:
        Dict with keys:
            - sequence_length: Total sequence length
            - enzymes_found: List of enzymes with cut sites
            - total_cut_sites: Total number of cut sites
            - sites: Dict of enzyme → positions
            - fragment_counts: Dict of enzyme → number of fragments produced
    """
    seq = dna_sequence.upper().replace('\n', '').replace('\r', '').replace(' ', '')
    sites = find_cut_sites(seq, enzymes)

    return {
        'sequence_length': len(seq),
        'enzymes_found': list(sites.keys()),
        'total_cut_sites': sum(len(v) for v in sites.values()),
        'sites': sites,
        'fragment_counts': {name: len(pos) + 1 for name, pos in sites.items()},
    }
