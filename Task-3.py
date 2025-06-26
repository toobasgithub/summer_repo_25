def read_fasta(filename):
    sequences = {}
# Error handling for input file
    try:
        with open(filename, 'r') as f:
            header = ''
            seq = ''
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header and seq:
                        sequences[header] = seq
                    header = line
                    seq = ''
                else:
                    if not header:
                        print("Error: Sequence found without header.")
                        return None
                    # To check sequence line only contain valid DNA characters.
                    VALID_CHARS = "ATCGN"
                    if not all(char.upper() in VALID_CHARS for char in line):
                        print(f"Error: Invalid character found in sequence on line {line_num} for header: {header}")
                        print(f"Invalid line content: '{line}'")
                        return None
                    seq += line.upper()
            if header and seq:
                sequences[header] = seq
        return sequences
    except FileNotFoundError: # Error handling
        print(f"Error: File '{filename}' not found.")
        return None
    except IOError:
        print(f"Error: Could not read file '{filename}'.")
    return None
# Store sequences with headers in a dictionary and validate the FASTA file format
# This step is embedded within the read_fasta function
# Filter by length
def filter_sequences(sequences, min_len):
    filtered = {}
    for header, seq in sequences.items():
        if len(seq) >= min_len:
            filtered[header] = seq
    return filtered
# Write to a new file
# Error handling for writing
def write_fasta(filename, sequences):
    try:
        with open(filename, 'w') as f:
            for header, seq in sequences.items():
                f.write(f"{header}\n")
                for i in range(0, len(seq), 60):  # Wrap at 60 chars
                    f.write(f"{seq[i:i+60]}\n")
        print(f"Filtered sequences saved to '{filename}'.")
    except IOError:
        print(f"Error: Could not write to file '{filename}'.")
# Main function to execute the script
def main():
    input_file = input("Enter input FASTA filename: ").strip()
    sequences = read_fasta(input_file)
    if sequences is None:
        return
    try:
        min_len = int(input("Enter minimum sequence length: ").strip())
    except ValueError:
        print("Error: Please enter a valid integer for sequence length.")
        return
    total = len(sequences)
    filtered = filter_sequences(sequences, min_len)
    passed = len(filtered)
    output_file = input("Enter output FASTA filename: ").strip()
    write_fasta(output_file, filtered)
# Display summary statistics
    print("\n--- Summary Statistics ---")
    print(f"Total sequences read: {total}")
    print(f"Sequences passing filter (>= {min_len} bp): {passed}")
    print(f"Sequences filtered out: {total - passed}")

if __name__ == "__main__":
    main()