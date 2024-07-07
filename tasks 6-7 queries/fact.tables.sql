CREATE TABLE fact_molecule_similarities (
    similarity_key INT PRIMARY KEY AUTO_INCREMENT,
    source_molecule_key INT,
    target_molecule_key INT,
    tanimoto_similarity DOUBLE PRECISION,
    has_duplicates_of_last_largest_score BOOLEAN,
    FOREIGN KEY (source_molecule_key) REFERENCES dim_molecules(molecule_key),
    FOREIGN KEY (target_molecule_key) REFERENCES dim_molecules(molecule_key)
);

INSERT INTO fact_molecule_similarities (source_molecule_key, target_molecule_key, tanimoto_similarity, has_duplicates_of_last_largest_score)
SELECT 
    src.molecule_key AS source_molecule_key,
    tgt.molecule_key AS target_molecule_key,
    ms.tanimoto_similarity,
    ms.has_duplicates_of_last_largest_score
FROM molecule_similarities ms
INNER JOIN dim_molecules src ON ms.source_chembl_id = src.chembl_id
INNER JOIN dim_molecules tgt ON ms.target_chembl_id = tgt.chembl_id;