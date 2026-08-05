from typing import Dict, Optional, Set

# Module-level constants
VALID_BASES: Set[str] = {'A', 'T', 'C', 'G'}


# Function 1: DNA Validator
def validate_dna_sequence(dna_sequence: str) -> bool:
    """
    Validate if a DNA sequence contains only valid bases (A, T, C, G).

    Args:
        dna_sequence: String containing DNA sequence

    Returns:
        True if valid.

    Raises:
        ValueError: If the sequence is empty or contains invalid characters.
    """
    dna_sequence = dna_sequence.upper()

    if not dna_sequence:
        raise ValueError("DNA sequence cannot be empty.")

    # More efficient way to check for invalid characters
    if not set(dna_sequence).issubset(VALID_BASES):
        invalid_chars = set(dna_sequence) - VALID_BASES
        raise ValueError(f"Invalid character(s) found in sequence: {', '.join(invalid_chars)}")

    return True


# Function 2: DNA Statistics
def calculate_dna_stats(dna_sequence: str) -> Optional[Dict]:
    """
    Calculate statistics for a DNA sequence.

    Args:
        dna_sequence: String containing DNA sequence

    Returns:
        Dictionary with length, counts, percentages, and gc_content
    """
    dna_sequence = dna_sequence.upper()
    length = len(dna_sequence)

    if length == 0:
        print("Error: DNA sequence cannot be empty.")
        return None # Or raise ValueError as in the validator

    counts: Dict[str, int] = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in dna_sequence:
        # This check is robust if the sequence wasn't pre-validated
        if base in counts:
            counts[base] += 1
        else:
            print(f"Warning: Invalid base '{base}' ignored.")

    percentages = {}
    for base in counts:
        percentages[base] = round((counts[base] / length) * 100, 2)

    gc_content = round(((counts['G'] + counts['C']) / length) * 100, 2)

    return {
        'length': length,
        'counts': counts,
        'percentages': percentages,
        'gc_content': gc_content
    }


# Function 3: FASTA File Reader
def read_fasta_file(filename: str) -> Optional[Dict[str, str]]:
    """
    Read a FASTA file and return sequences.

    Args:
        filename: Path to FASTA file

    Returns:
        Dictionary with sequence_name: sequence_data
    """
    sequences: Dict[str, str] = {}
    current_name: Optional[str] = None
    current_sequence = []

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()

                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = ''.join(current_sequence)
                    current_name = line[1:]
                    current_sequence = []
                elif line:
                    current_sequence.append(line.upper())
            if current_name:
                sequences[current_name] = ''.join(current_sequence)
        return sequences
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    

# Function 4: Transcribe DNA to RNA
def transcribe_dna_to_rna(dna_sequence: str) -> str:
    """
    Transcribes a DNA sequence into an RNA sequence by replacing 'T' with 'U'.

    Args:
        dna_sequence: A string containing the DNA sequence.

    Returns:
        A string representing the corresponding RNA sequence.
    """
    # Ensure sequence is valid before transcribing for robustness
    validate_dna_sequence(dna_sequence)
    return dna_sequence.upper().replace('T', 'U')
    