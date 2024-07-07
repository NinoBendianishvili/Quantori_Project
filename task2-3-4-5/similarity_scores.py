import os
import io
import csv
import boto3
import json
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
S3_BUCKET_NAME = "de-school-2024-aws"
S3_INPUT_PREFIX = "final_task/input_files/"
S3_CHEMBL_PREFIX = "final_task/bendianishvili_nino/chembl-fingerprints/"
S3_OUTPUT_PREFIX = "final_task/bendianishvili_nino/similarity_scores/"
AWS_REGION = "us-east-2"
FINGERPRINT_RADIUS = 2
FINGERPRINT_NBITS = 2048

# AWS Credentials
AWS_ACCESS_KEY_ID = 'hiding not to violate'
AWS_SECRET_ACCESS_KEY = 'hiding not to violate'
AWS_SESSION_TOKEN = 'hiding not to violate'


# Initialize S3 client
session = boto3.Session(
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
    aws_session_token=AWS_SESSION_TOKEN,
    region_name=AWS_REGION
)
s3_client = session.client('s3')

def read_s3_file(bucket, key):
    obj = s3_client.get_object(Bucket=bucket, Key=key)
    print(f"Reading file: {key}")
    buffer = io.BytesIO(obj['Body'].read())
    if buffer.getbuffer().nbytes == 0:
        print(f"Warning: The file {key} is empty.")
        return pd.DataFrame()
    if key.endswith('.csv'):
        try:
            return pd.read_csv(buffer)
        except pd.errors.ParserError:
            buffer.seek(0)
            return pd.read_csv(buffer, error_bad_lines=False, warn_bad_lines=True, quoting=csv.QUOTE_NONE)
    elif key.endswith('.parquet'):
        return pd.read_parquet(buffer)
    else:
        raise ValueError(f"Unsupported file format: {key}")

def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, FINGERPRINT_RADIUS, nBits=FINGERPRINT_NBITS)

def fingerprint_to_numpy(fp, n_bits=FINGERPRINT_NBITS):
    if isinstance(fp, str):
        return np.array([int(bit) for bit in fp], dtype=np.uint8)
    return np.array(fp, dtype=np.uint8)

def compute_tanimoto_similarity(input_fp, chembl_array):
    input_fp = input_fp.reshape(1, -1)
    intersection = np.dot(input_fp, chembl_array.T)
    input_sum = input_fp.sum()
    chembl_sum = chembl_array.sum(axis=1)
    union = input_sum + chembl_sum - intersection
    return (intersection / union).flatten()

def process_input_molecule(input_id, input_fp, chembl_array, chembl_ids):
    similarities = compute_tanimoto_similarity(input_fp, chembl_array)
    scores = list(zip(chembl_ids, similarities))
    scores.sort(key=lambda x: x[1], reverse=True)  # Sort by similarity score in descending order
    return input_id, scores

def save_results_to_s3(input_id, scores, bucket, prefix):
    csv_buffer = io.StringIO()
    writer = csv.writer(csv_buffer)
    writer.writerow(['ChEMBL_ID', 'Similarity_Score'])
    writer.writerows(scores)
    
    s3_client.put_object(
        Bucket=bucket,
        Key=f"{prefix}{input_id}_similarity_scores.csv",
        Body=csv_buffer.getvalue()
    )
    print(f"Saved results for input molecule: {input_id}")

def load_input_fingerprints():
    input_fingerprints = []
    response = s3_client.list_objects_v2(Bucket=S3_BUCKET_NAME, Prefix=S3_INPUT_PREFIX)
    if 'Contents' not in response:
        print(f"No input files found with prefix: {S3_INPUT_PREFIX}")
        return input_fingerprints

    for obj in response['Contents']:
        try:
            chunk = read_s3_file(S3_BUCKET_NAME, obj['Key'])
            if chunk.empty:
                continue
            for _, row in chunk.iterrows():
                if 'SMILES' not in row or 'Molecule Name' not in row:
                    print(f"Warning: Missing 'SMILES' or 'Molecule Name' in row: {row}")
                    continue
                fingerprint = smiles_to_fingerprint(row['SMILES'])
                if fingerprint is not None:
                    fp_string = fingerprint.ToBitString()
                    input_fingerprints.append((row['Molecule Name'], fp_string))
                else:
                    print(f"Warning: Could not generate fingerprint for SMILES: {row['SMILES']}")
        except Exception as e:
            print(f"Error processing file {obj['Key']}: {str(e)}")

    print(f"Successfully loaded {len(input_fingerprints)} input fingerprints")
    return input_fingerprints

def load_chembl_fingerprints():
    all_chembl_fps = []
    all_chembl_ids = []
    
    response = s3_client.list_objects_v2(Bucket=S3_BUCKET_NAME, Prefix=S3_CHEMBL_PREFIX)
    if 'Contents' not in response:
        print(f"No ChEMBL files found with prefix: {S3_CHEMBL_PREFIX}")
        return np.array([]), []

    for obj in response['Contents']:
        chembl_chunk = read_s3_file(S3_BUCKET_NAME, obj['Key'])
        all_chembl_fps.extend(chembl_chunk['fingerprint'].tolist())
        all_chembl_ids.extend(chembl_chunk['chembl_id'].tolist())
    
    return np.vstack([fingerprint_to_numpy(fp) for fp in all_chembl_fps]), all_chembl_ids

def main():
    print("Loading input fingerprints...")
    input_fingerprints = load_input_fingerprints()
    
    print("Loading ChEMBL fingerprints...")
    chembl_array, chembl_ids = load_chembl_fingerprints()
    print(f"Loaded {len(chembl_ids)} ChEMBL fingerprints")

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for input_id, input_fp in input_fingerprints:
            future = executor.submit(process_input_molecule, input_id, fingerprint_to_numpy(input_fp), chembl_array, chembl_ids)
            futures.append(future)
        
        for future in as_completed(futures):
            try:
                input_id, scores = future.result()
                save_results_to_s3(input_id, scores, S3_BUCKET_NAME, S3_OUTPUT_PREFIX)
            except Exception as e:
                print(f"Error processing molecule: {str(e)}")

if __name__ == "__main__":
    main()