from extract_promoters import parse_gtf, get_genome_dict, extract_promoter_sequences
from extract_features import extract_features
from model import train_model
import pandas as pd
import numpy as np

# Step 1: Load GTF and Genome
gtf_file = "data/annotation.gtf"
fasta_file = "data/genome.fa"
expr_file = "data/expression_data.txt"

promoters = parse_gtf(gtf_file)
genome = get_genome_dict(fasta_file)
sequences = extract_promoter_sequences(promoters, genome)

# Step 2: Extract Features
features = extract_features(sequences)
features_df = pd.DataFrame(features)

# Step 3: Load Expression and Merge
def load_expression_data(expression_file):
    df = pd.read_csv(expression_file, sep='\t', comment='#', skiprows=2)
    df = df.rename(columns={"Description": "gene"})
    tissue_columns = df.columns[2:] 
    df["expression"] = df[tissue_columns].mean(axis=1)
    df["expression"] = np.log1p(df["expression"])
    return df[["gene", "expression"]]

expression_df = load_expression_data(expr_file)
merged = pd.merge(features_df, expression_df, on="gene")
# Step 4: Train Model
train_model(merged)
print("Model trained successfully. Pipeline complete.")

# Show result with predictions
print("Sample of final output with predictions:")
print(pd.read_csv("predicted_gene_activity.csv").head())
