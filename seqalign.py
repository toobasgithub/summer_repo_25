from Bio import pairwise2
from Bio.Seq import Seq
import sys

# Function to perform global alignment
def align_seq(seq1, seq2, match=2, mismatch=-2, gap_open=-2, gap_extend=-1):
    alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    best_alignment = alignments[0]

    print("The alignment score is: ", best_alignment.score)
    print("The alignment sequence 1 is: ", best_alignment.seqA)
    print("The alignment sequence 2 is: ", best_alignment.seqB)
    print("The start of alignment is: ", best_alignment.start)
    print("The end of alignment is: ", best_alignment.end)

    return best_alignment

# Function to calculate seq similarity
def similarity(alignment):
    seq1 = alignment.seqA
    seq2 = alignment.seqB
    start = alignment.start
    end = alignment.end
    aligned1 = seq1[start:end]
    aligned2 = seq2[start:end]

    matches = 0
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i] and aligned1[i] != '-':
            matches += 1

    length = end - start
    similarity = (matches / length) * 100 if length > 0 else 0
    print("Matches: ", matches)
    print("The similarity of the alignment is:", round(similarity, 2), '%')
    return similarity

# Function to calculate gap frequency
def gap_frequency(alignment):
    seq1 = alignment.seqA
    seq2 = alignment.seqB

    gaps_seq1 = seq1.count('-')
    gaps_seq2 = seq2.count('-')
    total_gaps = gaps_seq1 + gaps_seq2
    alignment_length = len(seq1)

    gap_freq = (total_gaps / alignment_length) * 100 if alignment_length > 0 else 0

    print("Gaps in sequence 1:", gaps_seq1)
    print("Gaps in sequence 2:", gaps_seq2)
    print("Total gaps:", total_gaps)
    print("Gap frequency:", round(gap_freq, 2), '%')
    return gap_freq

# Function to calculate conserved regions
def find_conserved_regions(alignment, threshold=20):
    seq1, seq2 = alignment.seqA, alignment.seqB
    conserved_regions = []
    start = None

    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != '-':
            if start is None:
                start = i
        else:
            if start is not None and i - start >= threshold:
                region = seq1[start:i]
                conserved_regions.append(region)
            start = None

    # Final region check 
    if start is not None and len(seq1) - start >= threshold:
        region = seq1[start:]
        conserved_regions.append(region)
    
    print(f"Conserved Regions in Seq1 (≥ {threshold} bp)")
    if conserved_regions:
        for idx, region in enumerate(conserved_regions, 1):
            print(f"Conserved Region {idx} in Seq1:")
            print(f"Sequence: {region}\n")
    else:
        print(f"No conserved region of exactly {threshold} bp found.")

    return conserved_regions

if __name__ == '__main__':
    seq1 = Seq('ATGCCGTATGGCCAATGCTAGCCGTAAATGGCCGGTTAAAATGCTACGTAGTCGTACGTACGCCGTAAATTGGCCCCGGTTAATTGCATGACGTACAATCGTAGCGATCGATCGACGTACGTGCGTAAACGTACGTAACGGTAAACGTAACG')
    seq2 = Seq('ATGCCGTATGGCCAATGCTAAGGGTCCCTAGAGTCCGATAGCCGTAAATTGGCCCCGGTTAATTGCATGACGTACTTAGGTGCCATCGAGATCCCGTACGTGCGTAAACGTACGTAACGGTAAACGTAACG')
    alignment_results = align_seq(seq1, seq2)
    similarity(alignment_results)
    gap_frequency(alignment_results)
    find_conserved_regions(alignment_results, threshold=20)
    