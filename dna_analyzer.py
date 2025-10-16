# Function 1: DNA Validator 
def validate_dna_sequence(dna_sequence):

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

# Test both functions
print("=== Testing Validator ===")
print(validate_dna_sequence("ATCG"))

print("\n=== Testing Statistics ===")
result1 = calculate_dna_stats("ATCGATCG")
print(result1)

# Additional tests
print("\n=== Test with unequal distribution ===")
result2 = calculate_dna_stats("AAAGGGCCCTT")
print(result2)

print("\n=== Test with very high GC content ===")
result3 = calculate_dna_stats("GCGCGCGC")
print(result3)

print("\n=== Test with invalid character ===")
result4 = calculate_dna_stats("ATCGXYZ")
print(result4)
