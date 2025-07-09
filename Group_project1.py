from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Restriction import EcoRI, Analysis
import pandas as pd
import numpy as np
from collections import Counter
import sys

# Store all results here
results = []
file_path = 'saccharomyces_cerevisiae.fasta'

# Load sequence
try:
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            seq_id = record.id
            seq = record.seq.upper()

            # Validate
            valid_nucleotides = {'A', 'T', 'G', 'C'}
            if not set(seq).issubset(valid_nucleotides):
                print(f"Warning: Sequence {seq_id} contains invalid nucleotides. Skipping.")
                continue  

            # GC Content
            try:
                gc = round(gc_fraction(seq) * 100, 2)
            except:
                gc = 0.0

            # EcoRI site count
            try:
                analysis = Analysis([EcoRI], seq)
                ecori_sites = len(analysis.full().get('EcoRI', []))
            except:
                ecori_sites = 0

            # Transcription and Translation
            try:
                # trim the sequence to a length that is a multiple of 3
                trimmed_len = len(seq) - (len(seq) % 3)
                coding_dna = seq[:trimmed_len]
                # Use the new 'coding_dna' variable for BOTH steps
                rna = coding_dna.transcribe()
                protein = coding_dna.translate(to_stop=True)
            except Exception:
                        rna = "ERROR"
                        protein = "ERROR"  
            
            # Codon usage
            codon_freq = Counter()
            try:
                coding_seq = seq[:len(seq) - len(seq) % 3]
                for i in range(0, len(coding_seq), 3):
                    codon = str(coding_seq[i:i+3])
                    codon_freq[codon] += 1
            except:
                codon_freq = {}

            # Store everything
            results.append({
                'Sequence_ID': seq_id,
                'Length': len(seq),
                'GC_Content': gc,
                'EcoRI_Sites': ecori_sites,
                'RNA': str(rna)[:50],
                'Protein': str(protein)[:50],
                'Codon_Usage': dict(codon_freq)
            })

except FileNotFoundError:
    print("FASTA file not found.")
except Exception as e:
    print(f"Error while reading sequences: {e}")

try:
    # Create DataFrame
    df = pd.DataFrame(results)
    if df.empty:
        print("No valid sequences found in the file. Cannot proceed with analysis.")
        sys.exit()
    # Save to CSV
    df.to_csv('sequence_analysis.csv', index=False)
    print("Results saved to sequence_analysis.csv successfully.")
except Exception as e:
    print("Error writing to CSV:", e)
   
# NumPy: Sequence length distribution 
lengths = df['Length'].values
mean_length = np.mean(lengths)
std_dev = np.std(lengths)
print(f"Mean Sequence Length: {mean_length:.2f}")
print(f"Standard Deviation: {std_dev:.2f}")
    
# Bar Plot: GC Content per sequence
import matplotlib.pyplot as plt
def plot_gc_content_bar(df):
    plt.figure(figsize=(10, 6))
    plt.bar(df['Sequence_ID'], df['GC_Content'], color='skyblue')
    plt.xlabel('Sequence ID')
    plt.ylabel('GC Content (%)')
    plt.title('GC Content per Sequence')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('gc_content_bar.png')
    plt.close()
print("Saved: gc_content_bar.png")
plot_gc_content_bar(df)

# heatmap: codon usuage frequencies
import seaborn as sns
def plot_codon_heatmap(df):
    codon_usage = df['Codon_Usage'].iloc[0]
    usage_df = pd.DataFrame.from_dict(codon_usage, orient='index', columns=['Frequency'])
    usage_df = usage_df.reset_index().rename(columns={'index': 'Codon'})
    usage_df = usage_df.pivot_table(index='Codon', values='Frequency')
    plt.figure(figsize=(10, 8))
    sns.heatmap(usage_df, annot=True, fmt=".2f", cmap='YlGnBu', cbar=True)
    plt.title('Codon Usage Heatmap (First Sequence)')
    plt.savefig('codon_heatmap.png')
    plt.close()
    print("Saved: codon_heatmap.png")
plot_codon_heatmap(df)

# histogram: matplotlib distribution of sequence length    
def plot_length_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df['Length'], bins=10, color='orange', edgecolor='black')
    plt.xlabel('Sequence Length')
    plt.ylabel('Count')
    plt.title('Distribution of Sequence Lengths')
    plt.savefig('length_histogram.png')
    plt.close()
print("Saved: length_histogram.png")
plot_length_histogram(df)

# Plotly Scatter Plot: Plot GC content vs. EcoRI sites
import plotly.express as px
def plot_gc_vs_ecori(df):
    fig = px.scatter(
        df,
        x='GC_Content',
        y='EcoRI_Sites',
        size='Length',
        color='Length',
        hover_name='Sequence_ID',
        title='GC Content vs EcoRI Sites',
        labels={'GC_Content': 'GC Content (%)', 'EcoRI_Sites': 'EcoRI Sites'}
    )
    fig.write_html('restriction_scatter.html')
    print("Saved: restriction_scatter.html")
plot_gc_vs_ecori(df)









