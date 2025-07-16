import gzip
import csv
import re

def parse_vcf_line(line):
    fields = line.strip().split('\t')
    chrom, pos, vid, ref, alt, qual, filter_, info = fields[:8]
    
    clnsig_match = re.search(r'CLNSIG=([^;]+)', info)
    clnsig = clnsig_match.group(1) if clnsig_match else None
    
    af_match = re.search(r'AF_ESP=([^;]+)', info)
    allele_freq = float(af_match.group(1)) if af_match else 0.0
    
    return chrom, pos, ref, alt, clnsig, allele_freq

def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT'

def label_from_clnsig(clnsig):
    if clnsig is None: return None
    if 'Pathogenic' in clnsig or 'Likely_pathogenic' in clnsig: return 'deleterious'
    if 'Benign' in clnsig or 'Likely_benign' in clnsig: return 'neutral'
    return None

def one_hot_base(base):
    return [int(base == b) for b in 'ACGT']

def extract_features(input_vcf_gz, output_csv):
    with gzip.open(input_vcf_gz, 'rt') as vcf, open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['chrom', 'pos', 'ref_A', 'ref_C', 'ref_G', 'ref_T', 'alt_A', 'alt_C', 'alt_G', 'alt_T', 'allele_freq', 'label'])
        
        for line in vcf:
            if line.startswith('#'): continue
                
            chrom, pos, ref, alt, clnsig, allele_freq = parse_vcf_line(line)
            
            if not is_snv(ref, alt): continue
            label = label_from_clnsig(clnsig)
            if label is None: continue
                
            ref_oh = one_hot_base(ref)
            alt_oh = one_hot_base(alt)
            
            writer.writerow([chrom, pos] + ref_oh + alt_oh + [allele_freq] + [label])
            
    print(f"Features extracted to '{output_csv}'")