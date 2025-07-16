import os
import requests

def download_clinvar(url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz", output_file="clinvar.vcf.gz"):
    """
    Downloads the ClinVar VCF file from the given URL and saves it locally.
    """
    if os.path.exists(output_file):
        print(f"'{output_file}' already exists. Skipping download.")
        return
    try:
        print(f"Downloading {url} ...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print(f"Downloaded and saved as '{output_file}'")
    except Exception as e:
        print(f"Error downloading file: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)