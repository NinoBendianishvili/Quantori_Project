CREATE TABLE dim_molecules (
    molecule_key INT PRIMARY KEY AUTO_INCREMENT,
    chembl_id VARCHAR(50) UNIQUE NOT NULL,
    molecule_type VARCHAR(50),
    mw_freebase VARCHAR(50),
    alogp VARCHAR(50),
    psa VARCHAR(50),
    cx_logp VARCHAR(50),
    molecular_species VARCHAR(50),
    full_mwt VARCHAR(50),
    aromatic_rings DOUBLE PRECISION,
    heavy_atoms DOUBLE PRECISION
);


INSERT INTO dim_molecules (chembl_id, molecule_type, mw_freebase, alogp, psa, cx_logp, molecular_species, full_mwt, aromatic_rings, heavy_atoms)
SELECT DISTINCT 
    cp.chembl_id, 
    md.molecule_type, 
    cp.mw_freebase, 
    cp.alogp, 
    cp.psa, 
    cp.cx_logd, 
    cp.molecular_species, 
    cp.full_mwt, 
    cp.aromatic_rings, 
    cp.heavy_atoms
FROM staging_compound_properties cp
INNER JOIN staging_molecule_dictionary md ON cp.chembl_id = md.chembl_id
INNER JOIN (
    SELECT source_chembl_id AS chembl_id FROM molecule_similarities
    UNION
    SELECT target_chembl_id FROM molecule_similarities
) s ON cp.chembl_id = s.chembl_id;