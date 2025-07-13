from Bio import SeqIO
import pandas as pd

def parse_gtf(gtf_file):
    promoters = []
    with open(gtf_file) as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "transcript":
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                info = fields[8]
                gene_name_match = [x for x in info.split(';') if 'gene_name' in x]
                if not gene_name_match:
                    continue  # skip if no gene name
                gene_name = gene_name_match[0].split('"')[1]

                if strand == '+':
                    promoter_start = max(1, start - 1000)
                    promoter_end = start
                else:
                    promoter_start = end
                    promoter_end = end + 1000

                promoters.append((chrom, promoter_start, promoter_end, strand, gene_name))
    return promoters

def get_genome_dict(fasta_path):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    return genome

def extract_promoter_sequences(promoter_regions, genome):
    sequences = []
    for chrom, start, end, strand, gene in promoter_regions:
        if chrom in genome:
            seq = genome[chrom].seq[start:end]
            if strand == '-':
                seq = seq.reverse_complement()
            sequences.append((gene, str(seq)))
    return sequences