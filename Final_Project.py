import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import gzip
import re
import os
import numpy as np 

def main():
    # --- Step 1 & 2: Choose Organism and Curate Data ---
    GFF_FILE = 'Arabidopsis_thaliana.TAIR10.61.gff3.gz'
    VCF_FILE = 'arabidopsis_thaliana.vcf.gz'
    EXPRESSION_FILE = 'gene_expression.tsv'
    ORGANISM_NAME = 'Arabidopsis thaliana'
    OUTPUT_DIR = 'output'
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # --- Step 3: Parse Genome Annotation (GFF file) ---
    print("Parsing GFF file for exonic regions...")
    gene_exons = {}
    with gzip.open(GFF_FILE, 'rt', errors='ignore') as gff_handle:
        for line in gff_handle:
            if line.startswith('#'):
                continue 
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  
            if parts[2] == 'exon':
                match = re.search(r'Parent=transcript:(AT\dG\d{5})', parts[8])
                if match:
                    gene_id = match.group(1)
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                if gene_id not in gene_exons:
                    gene_exons[gene_id] = {'chrom': chrom, 'exons': []}
                gene_exons[gene_id]['exons'].append((start, end))
    print(f"Parsed {len(gene_exons)} genes with exon coordinates.")
    
    # --- Step 4: Parse SNPs from VCF file ---
    print("Loading SNPs from VCF file...")
    snps = pd.read_csv(VCF_FILE, sep='\t', comment='#', header=None, usecols=[0, 1], names=['chrom', 'pos'], dtype={'chrom': str})
    print(f"Loaded {len(snps)} SNPs.")
    
    # --- Step 5: Map SNPs to Exons & Step 6: Calculate SNP Density ---
    print("Mapping SNPs to exons and calculating SNP density...")
    print("Optimizing by pre-grouping SNPs by chromosome...")
    snps_by_chrom = {chrom: group['pos'].to_numpy() for chrom, group in snps.groupby('chrom')}
    results = []
    gene_count = 0
    total_genes = len(gene_exons)
    for gene_id, data in gene_exons.items():
        gene_count += 1
        if gene_count % 1000 == 0:
            print(f"   -> Processing gene {gene_count}/{total_genes}...")
        chrom = data['chrom']
        exons = data['exons']
        total_exon_length = sum(end - start for start, end in exons)
        if total_exon_length == 0:
            continue
        snp_count = 0
        chrom_snps_pos = snps_by_chrom.get(chrom)
        if chrom_snps_pos is not None:
            for start, end in exons:
                count = np.count_nonzero((chrom_snps_pos >= start) & (chrom_snps_pos <= end))
                snp_count += count 
        density = snp_count / total_exon_length
        results.append({'gene_id': gene_id, 'snp_density': density})
    snp_density_df = pd.DataFrame(results)
    print(f"Calculated SNP density for {len(snp_density_df)} genes.")

    # --- Step 7: Correlate with Gene Expression ---
    print("Merging SNP density with gene expression data...")
    expression_data = pd.read_csv(EXPRESSION_FILE, sep='\t')
    expression_data = expression_data.rename(columns={'GeneID': 'gene_id', 'ExpressionLevel': 'expression_level'})
    merged_data = pd.merge(snp_density_df, expression_data, on='gene_id', how='inner')
    merged_data = merged_data.dropna()
    merged_data['log_expression'] = np.log10(merged_data['expression_level'] + 1)
    unmatched = set(expression_data['gene_id']) - set(merged_data['gene_id'])
    print(f"Unmatched gene IDs (not found in SNP data): {list(unmatched)[:10]}")
    print(f"Matched: {len(merged_data)}, Unmatched: {len(unmatched)}")
    print(f"Final merged dataset contains {len(merged_data)} genes.")

    # --- Step 8: Create Figure to Visualize Relationship ---
    print("Creating scatter plot...")
    plt.figure(figsize=(10, 7))
    sns.scatterplot(x='log_expression', y='snp_density', data=merged_data, alpha=0.4, s=15)
    plt.xlabel('Gene Expression Level (log10 scale)')
    plt.ylabel('SNP Density in Exons (SNPs per base pair)')
    plt.title(f'SNP Density vs. Gene Expression in {ORGANISM_NAME}')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    output_plot_path = os.path.join(OUTPUT_DIR, 'snp_density_vs_expression.png')
    plt.savefig(output_plot_path)
    print(f"Plot saved to {output_plot_path}")
    # plt.show()

    # --- Step 9: Perform Pearson Correlation Test ---
    print("Performing Pearson correlation test...")
    if len(merged_data) > 1:
        r, p_value = pearsonr(merged_data['log_expression'], merged_data['snp_density'])
        print(f"Pearson Correlation between log(Expression) and SNP Density:")
        print(f"Correlation coefficient (r) = {r:.4f}")
        print(f"P-value = {p_value:.4e}")
        if p_value < 0.05:
            print("Result: The correlation is statistically significant (p < 0.05).")
        else:
            print("Result: The correlation is not statistically significant (p >= 0.05).")
    else:
        print("Not enough data to perform correlation test.")

    # --- Save final data and complete ---
    output_csv_path = os.path.join(OUTPUT_DIR, 'final_analysis_data.csv')
    merged_data.to_csv(output_csv_path, index=False)
    print(f"Analysis complete. Final data saved to {output_csv_path}")

if __name__ == "__main__":
    main()