import psycopg2
import requests
import json
import time
from multiprocessing import Pool, cpu_count
def fetch_data_and_column_names(endpoint, limit=1000, offset=0, max_retries=5):
    base_url = 'https://www.ebi.ac.uk/chembl/api/data'
    url = f"{base_url}/{endpoint}.json"
    params = {'limit': limit, 'offset': offset}
    retries = 0

    while retries < max_retries:
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            key = endpoint + 's'
            if key in data:
                columns = list(data[key][0].keys())
                print(f"Fetched {len(data[key])} records from offset {offset}")
                return data[key], columns
            else:
                raise KeyError(f"The expected key '{key}' is not present in the response.")
        except (requests.exceptions.RequestException, KeyError) as e:
            retries += 1
            print(f"Error fetching data (attempt {retries}/{max_retries}) from offset {offset}: {e}")
            time.sleep(2 ** retries)
    print(f"Skipping offset {offset} after {max_retries} retries")
    return [], []

def create_table(conn_params, table_name, columns):
    try:
        with psycopg2.connect(**conn_params) as conn:
            with conn.cursor() as cur:
                columns_str = ", ".join(f"{col} TEXT" for col in columns)
                create_table_sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({columns_str});"
                cur.execute(create_table_sql)
                conn.commit()
    except psycopg2.Error as e:
        print(f"Error creating table {table_name}: {e}")

def serialize_data(record):
    return {k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in record.items()}
def insert_into_postgres_chunk(conn_params, table_name, chunk):
    try:
        with psycopg2.connect(**conn_params) as conn:
            with conn.cursor() as cur:
                for record in chunk:
                    serialized_record = serialize_data(record)
                    columns_str = ", ".join(serialized_record.keys())
                    placeholders = ", ".join(["%s"] * len(serialized_record))
                    insert_sql = f"INSERT INTO {table_name} ({columns_str}) VALUES ({placeholders});"
                    cur.execute(insert_sql, tuple(serialized_record.values()))
                conn.commit()
    except psycopg2.Error as e:
        print(f"Error inserting into table {table_name}: {e}")
def fetch_worker(params):
    endpoint, limit, offset, max_retries = params
    data, columns = fetch_data_and_column_names(endpoint, limit, offset, max_retries)
    return data, columns, offset


def ingest_chembl(endpoint, max_processes=None, max_retries=5, batch_size=5, delay_between_batches=10):
    try:
        conn_params = {
            'host': "de-database.ccpquxzb7z8t.us-east-2.rds.amazonaws.com",
            'database': "postgres",
            'user': "nbendianishvili",
            'password': "c56f257b"
        }
        offset = 0
        limit = 1000
        initial_data, columns = fetch_data_and_column_names(endpoint, limit, offset, max_retries)
        if initial_data:
            create_table(conn_params, endpoint, columns)
        if max_processes is None:
            max_processes = cpu_count() 

        with Pool(processes=max_processes) as pool:
            while True:
                params_list = [(endpoint, limit, offset + i * limit, max_retries) for i in range(batch_size)]
                results = pool.map(fetch_worker, params_list)

                for data, columns, current_offset in results:
                    if data:
                        insert_into_postgres_chunk(conn_params, endpoint, data)
                if len(initial_data) < limit:
                    break
                offset += limit * batch_size
                time.sleep(delay_between_batches)
    except Exception as e:
        print("Error during ingestion process:", e)

if __name__ == "__main__":
    ingest_chembl('molecule')
    ingest_chembl("chembl_id_lookup")
