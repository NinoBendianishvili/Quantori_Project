import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import boto3
import pyarrow.parquet as pq
import pyarrow as pa
from sqlalchemy import create_engine, text
import multiprocessing as mp
from functools import partial
import pickle

# Configuration
DATABASE_CONFIG = {
    'dbname': 'postgres',
    'user': 'nbendianishvili',
    'password': 'c56f257b',
    'host': 'de-database.ccpquxzb7z8t.us-east-2.rds.amazonaws.com',
    'port': '5432'
}
S3_BUCKET_NAME = "de-school-2024-aws"
S3_FOLDER_NAME = "final_task/bendianishvili_nino/chembl-fingerprints/"
AWS_REGION = "us-east-2"
FINGERPRINT_RADIUS = 2
FINGERPRINT_NBITS = 2048
CHUNK_SIZE = 50000  
NUM_PROCESSES = mp.cpu_count() 
os.environ['AWS_ACCESS_KEY_ID'] = 'hiding not to violate'
os.environ['AWS_SECRET_ACCESS_KEY'] = 'hiding not to violate'
os.environ['AWS_DEFAULT_REGION'] = AWS_REGION   
# AWS S3 Client
s3_client = boto3.client('s3', region_name=AWS_REGION)

def fetch_compound_structures(chunk_size, engine):
    query = text("SELECT chembl_id, canonical_smiles FROM staging_compound_structures")
    with engine.connect() as conn:
        result = conn.execution_options(stream_results=True).execute(query)
        while True:
            chunk = result.fetchmany(chunk_size)
            if not chunk:
                break
            yield pd.DataFrame(chunk, columns=['chembl_id', 'canonical_smiles'])

def compute_morgan_fingerprint(row, radius=FINGERPRINT_RADIUS, nBits=FINGERPRINT_NBITS):
    smiles = row['canonical_smiles']
    chembl_id = row['chembl_id']
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        return {"chembl_id": chembl_id, "fingerprint": fingerprint.ToBitString()}
    else:
        return None

def process_chunk(chunk):

    with mp.Pool(NUM_PROCESSES) as pool:
        results = pool.map(compute_morgan_fingerprint, [row for _, row in chunk.iterrows()])
    return [r for r in results if r is not None]

def save_chunk_to_parquet(chunk, file_path):
    table = pa.Table.from_pandas(pd.DataFrame(chunk))
    pq.write_table(table, file_path)

def upload_to_s3(file_path, s3_bucket, s3_key):
    print("started uploading")
    try:
        s3_client = boto3.client(
            's3',
            aws_access_key_id='hiding not to violate',
            aws_secret_access_key='hiding not to violate',
            aws_session_token='hiding not to violate')
        s3_client.upload_file(file_path, s3_bucket, s3_key)
        print(f"Successfully uploaded {file_path} to {s3_bucket}/{s3_key}")
    except Exception as e:
        print(f"Failed to upload {file_path} to {s3_bucket}/{s3_key}: {e}")
        
def main():
    chunk_number = 0
    engine = create_engine(f"postgresql+psycopg2://{DATABASE_CONFIG['user']}:{DATABASE_CONFIG['password']}@{DATABASE_CONFIG['host']}:{DATABASE_CONFIG['port']}/{DATABASE_CONFIG['dbname']}")

    print("Fetching and processing compound structures in chunks")
    for chunk in fetch_compound_structures(CHUNK_SIZE, engine):
        print(f"Processing chunk {chunk_number}")
        
        all_fingerprints = process_chunk(chunk)

        fingerprints_file_path = f"morgan_fingerprints_chunk_{chunk_number}.parquet"
        save_chunk_to_parquet(all_fingerprints, fingerprints_file_path)

        s3_key = f"{S3_FOLDER_NAME}morgan_fingerprints_chunk_{chunk_number}.parquet"
        upload_to_s3(fingerprints_file_path, S3_BUCKET_NAME, s3_key)
        
        # Clean up the temporary file
        os.remove(fingerprints_file_path)

        chunk_number += 1

if __name__ == "__main__":
    main()