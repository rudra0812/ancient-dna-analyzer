# Function name: validate_dna_sequence
# Input: A string (DNA sequence)
# Output: True if valid, False if invalid

#Rules:
# Only A, T, C, G letters allowed (uppercase)
# - Should handle empty strings
#- Print helpful error message if invalid 

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

    # Test cases
print(validate_dna_sequence("ATCG"))        
print(validate_dna_sequence("ATXG"))        
print(validate_dna_sequence(""))            
print(validate_dna_sequence("atcg"))        
print(validate_dna_sequence("ATCG123"))
