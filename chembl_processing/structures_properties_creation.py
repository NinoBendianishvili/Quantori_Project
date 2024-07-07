from sqlalchemy import create_engine, MetaData, Table, Column, String, ForeignKey, inspect, Float, Integer, Boolean
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.exc import IntegrityError
import json

# Database connection configuration
DATABASE_URI = 'postgresql+psycopg2://nbendianishvili:c56f257b@de-database.ccpquxzb7z8t.us-east-2.rds.amazonaws.com/postgres'

# Create a new database engine instance
engine = create_engine(DATABASE_URI)

# Initialize metadata
metadata = MetaData()

# Reflect the existing molecule table
molecule_table = Table('molecule', metadata, autoload_with=engine)

# Define the new molecule_structures table
molecule_structures_table = Table(
    'molecule_structures', metadata,
    Column('molecule_chembl_id', String, ForeignKey('molecule.molecule_chembl_id'), primary_key=True),
    Column('canonical_smiles', String),
    Column('molfile', String),
    Column('standard_inchi', String),
    Column('standard_inchi_key', String)
)

# Define the new molecule_properties table
molecule_properties_table = Table(
    'molecule_properties', metadata,
    Column('molecule_chembl_id', String, ForeignKey('molecule.molecule_chembl_id'), primary_key=True),
    Column('alogp', Float),
    Column('aromatic_rings', Integer),
    Column('cx_logd', Float),
    Column('cx_logp', Float),
    Column('cx_most_apka', Float),
    Column('cx_most_bpka', Float),
    Column('full_molformula', String),
    Column('full_mwt', Float),
    Column('hba', Integer),
    Column('hba_lipinski', Integer),
    Column('hbd', Integer),
    Column('hbd_lipinski', Integer),
    Column('heavy_atoms', Integer),
    Column('molecular_species', String),
    Column('mw_freebase', Float),
    Column('mw_monoisotopic', Float),
    Column('np_likeness_score', Float),
    Column('num_lipinski_ro5_violations', Integer),
    Column('num_ro5_violations', Integer),
    Column('psa', Float),
    Column('qed_weighted', Float),
    Column('ro3_pass', String),
    Column('rtb', Integer)
)

# Create the tables if they don't exist
inspector = inspect(engine)
if not inspector.has_table('molecule_structures'):
    molecule_structures_table.create(engine)
if not inspector.has_table('molecule_properties'):
    molecule_properties_table.create(engine)

# Create a session
Session = sessionmaker(bind=engine)

def upsert_data(session, table, batch):
    insert_stmt = insert(table).values(batch)
    upsert_stmt = insert_stmt.on_conflict_do_update(
        index_elements=['molecule_chembl_id'],
        set_={c.name: insert_stmt.excluded[c.name] for c in table.columns if c.name != 'molecule_chembl_id'}
    )
    session.execute(upsert_stmt)

def extract_and_insert_data(batch_size=1000):
    session = Session()
    try:
        offset = 0
        while True:
            # Query the molecule table in batches
            query = session.query(molecule_table.c.molecule_chembl_id, 
                                  molecule_table.c.molecule_structures,
                                  molecule_table.c.molecule_properties)\
                .offset(offset).limit(batch_size)
            results = query.all()

            if not results:
                break

            structures_batch = []
            properties_batch = []
            for result in results:
                molecule_chembl_id = result.molecule_chembl_id
                molecule_structures = json.loads(result.molecule_structures) if result.molecule_structures else None
                molecule_properties = json.loads(result.molecule_properties) if result.molecule_properties else None

                if molecule_structures:
                    structures_batch.append({
                        'molecule_chembl_id': molecule_chembl_id,
                        'canonical_smiles': molecule_structures.get('canonical_smiles'),
                        'molfile': molecule_structures.get('molfile'),
                        'standard_inchi': molecule_structures.get('standard_inchi'),
                        'standard_inchi_key': molecule_structures.get('standard_inchi_key')
                    })

                if molecule_properties:
                    properties_batch.append({
                        'molecule_chembl_id': molecule_chembl_id,
                        'alogp': molecule_properties.get('alogp'),
                        'aromatic_rings': molecule_properties.get('aromatic_rings'),
                        'cx_logd': molecule_properties.get('cx_logd'),
                        'cx_logp': molecule_properties.get('cx_logp'),
                        'cx_most_apka': molecule_properties.get('cx_most_apka'),
                        'cx_most_bpka': molecule_properties.get('cx_most_bpka'),
                        'full_molformula': molecule_properties.get('full_molformula'),
                        'full_mwt': molecule_properties.get('full_mwt'),
                        'hba': molecule_properties.get('hba'),
                        'hba_lipinski': molecule_properties.get('hba_lipinski'),
                        'hbd': molecule_properties.get('hbd'),
                        'hbd_lipinski': molecule_properties.get('hbd_lipinski'),
                        'heavy_atoms': molecule_properties.get('heavy_atoms'),
                        'molecular_species': molecule_properties.get('molecular_species'),
                        'mw_freebase': molecule_properties.get('mw_freebase'),
                        'mw_monoisotopic': molecule_properties.get('mw_monoisotopic'),
                        'np_likeness_score': molecule_properties.get('np_likeness_score'),
                        'num_lipinski_ro5_violations': molecule_properties.get('num_lipinski_ro5_violations'),
                        'num_ro5_violations': molecule_properties.get('num_ro5_violations'),
                        'psa': molecule_properties.get('psa'),
                        'qed_weighted': molecule_properties.get('qed_weighted'),
                        'ro3_pass': molecule_properties.get('ro3_pass'),
                        'rtb': molecule_properties.get('rtb')
                    })

            if structures_batch:
                upsert_data(session, molecule_structures_table, structures_batch)
            if properties_batch:
                upsert_data(session, molecule_properties_table, properties_batch)

            session.commit()
            offset += batch_size
            print(f"Processed {offset} records")

    except Exception as e:
        print(f"An error occurred: {e}")
        session.rollback()
    finally:
        session.close()

if __name__ == "__main__":
    extract_and_insert_data()