import re
import numpy as np
from collections import Counter

def kmer_counts(seq, k=3):
    return Counter([seq[i:i+k] for i in range(len(seq)-k+1)])

def gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) if len(seq) > 0 else 0

def at_content(seq):
    a = seq.count('A')
    t = seq.count('T')
    return (a + t) / len(seq) if len(seq) > 0 else 0

def has_tata_box(seq):
    return 'TATAAA' in seq

def dinucleotide_repeat_score(seq):
    # Count repeated dinucleotide motifs (e.g., ATATAT or GCGCGC)
    count = 0
    for dinuc in ['AT', 'TA', 'CG', 'GC', 'AG', 'GA', 'CT', 'TC']:
        pattern = f"({dinuc}){{3,}}"  # 3 or more repeats
        if re.search(pattern, seq):
            count += 1
    return count

def extract_features(sequences):
    feature_list = []
    for gene, seq in sequences:
        features = {
            "gene": gene,
            "gc_content": gc_content(seq),
            "at_content": at_content(seq),
            "tata_box": int(has_tata_box(seq)),
            "dinuc_repeat_score": dinucleotide_repeat_score(seq),
        }

        # Add 3-mer frequencies
        kmer_freqs = kmer_counts(seq, 3)
        for kmer in kmer_freqs:
            features[f'kmer_{kmer}'] = kmer_freqs[kmer]

        feature_list.append(features)
    return feature_list