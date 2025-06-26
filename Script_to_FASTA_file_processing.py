import sys
import csv
#Step1.Function to calculate GC content
def calculate_gc_content(seq):
    if len(seq) == 0:
        return 0
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)
#Step2.Function to validate DNA sequence
def is_valid_dna(seq):
    return all(base in 'ATGC' for base in seq)
#Step3.Main script logic
if __name__ == "__main__":
    # Input file check
    if len(sys.argv) != 2:
        print("Usage: python Script_to_FASTA_file_processing.py <input_fasta>")
        sys.exit(1)
    input_file = sys.argv[1]
#Step4.Read FASTA file and store sequences in a dictionary
sequences = {}
try:
    with open(input_file, 'r') as f:
        current_id = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ''
            else:
                sequences[current_id] += line.upper()
except FileNotFoundError:
    print(f"Error: File '{input_file}' not found.")
    sys.exit(1)
#Step5.Use set to find unique nucleotides
unique_nucleotides = set()
for seq in sequences.values():
    unique_nucleotides.update(seq)
print("Unique nucleotides in all sequences:", sorted(unique_nucleotides))
#Step6.Analyze each sequence and prepare results
results = []
for seq_id, seq in sequences.items():
    gc_content = calculate_gc_content(seq)
    valid = is_valid_dna(seq)
    result = {
        'ID': seq_id,
        'Length': len(seq),
        'GC_Content': round(gc_content, 2),
        'Is_Valid': valid
        }
    results.append(result)
#Step7.Save results to CSV file
output_file = "fasta_analysis_output.csv"
try:
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['ID', 'Length', 'GC_Content', 'Is_Valid']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"\nAnalysis complete. Results saved to '{output_file}'")
except IOError:
    print(f"Error: Could not write to '{output_file}'.")