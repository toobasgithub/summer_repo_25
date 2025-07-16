import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

from data_download import download_clinvar
from feature_extraction import extract_features
from train import run_training_pipeline

def main():
    """
    This is the main controller that runs the entire pipeline.
    """
    # Define file paths
    vcf_gz_file = 'clinvar.vcf.gz'
    features_csv_file = 'snv_features.csv'

    # Download Data 
    print("DOWNLOADING DATA")
    download_clinvar(output_file=vcf_gz_file)

    # Extract Features 
    print("EXTRACTING FEATURES")
    extract_features(input_vcf_gz=vcf_gz_file, output_csv=features_csv_file)

    # Train and Evaluate Model
    print("TRAINING AND EVALUATING MODEL")
    run_training_pipeline(csv_file=features_csv_file)
    
    print("PIPELINE COMPLETED SUCCESSFULLY!")

if __name__ == "__main__":
    main()