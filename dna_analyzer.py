# Function 1: DNA Validator 
def validate_dna_sequence(dna_sequence):
    """
    Validate if a DNA sequence contains only valid bases (A, T, C, G).
    
    Args:
        dna_sequence: String containing DNA sequence
        
    Returns:
        True if valid, False otherwise
    """
    dna_sequence = dna_sequence.upper()
    valid_bases = {'A', 'T', 'C', 'G'}

    if not dna_sequence:
        print("Error: DNA sequence cannot be empty.")
        return False

    for base in dna_sequence:
        if base not in valid_bases:
            print(f"Error: Invalid character '{base}' found in the sequence.")
            return False
        
    print("DNA sequence is valid.")
    return True


# Function 2: DNA Statistics 
def calculate_dna_stats(dna_sequence):
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
        return None
    
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in dna_sequence:
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
def read_fasta_file(filename):
    """
    Read a FASTA file and return sequences.
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Dictionary with sequence_name: sequence_data
    """
    sequences = {}
    current_name = None
    current_sequence = []
    
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_name:
                        sequences[current_name] = ''.join(current_sequence)
                    
                    # Start new sequence
                    current_name = line[1:]  # Remove '>'
                    current_sequence = []
                    
                elif line:
                    # Sequence line
                    current_sequence.append(line.upper())
            
            # Save last sequence
            if current_name:
                sequences[current_name] = ''.join(current_sequence)
                
        return sequences
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None


# =============================================================================
# TESTS
# =============================================================================

if __name__ == "__main__":
    # Test 1: Validator
    print("=== Testing Validator ===")
    print(validate_dna_sequence("ATCG"))
    
    # Test 2: Basic Statistics
    print("\n=== Testing Statistics ===")
    result1 = calculate_dna_stats("ATCGATCG")
    print(result1)
    
    # Test 3: Unequal distribution
    print("\n=== Test with unequal distribution ===")
    result2 = calculate_dna_stats("AAAGGGCCCTT")
    print(result2)
    
    # Test 4: High GC content
    print("\n=== Test with very high GC content ===")
    result3 = calculate_dna_stats("GCGCGCGC")
    print(result3)
    
    # Test 5: Invalid characters
    print("\n=== Test with invalid character ===")
    result4 = calculate_dna_stats("ATCGXYZ")
    print(result4)
    
    # Test 6: FASTA reader
    print("\n=== Testing FASTA Reader ===")
    sequences = read_fasta_file("sample_sequences.fasta")
    
    if sequences:
        print(f"Found {len(sequences)} sequences:\n")
        
        for name, seq in sequences.items():
            print(f"Sequence: {name}")
            print(f"Length: {len(seq)}")
            
            # Analyze each sequence
            stats = calculate_dna_stats(seq)
            if stats:
                print(f"GC Content: {stats['gc_content']}%")
            print("-" * 50)