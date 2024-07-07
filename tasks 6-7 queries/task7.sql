CREATE VIEW mol_similarity_averages AS (
    SELECT 
        source_chembl_id AS origin_molecule_id,
        ROUND(AVG(tanimoto_similarity)::numeric, 4) AS mean_similarity
    FROM 
        molecule_similarities
    GROUP BY 
        source_chembl_id
);

-- B
CREATE VIEW molecule_lipophilicity_comparisons AS
WITH molecule_characteristics AS (
    SELECT 
        chembl_id AS molecule_identifier, 
        alogp AS lipophilicity_index
    FROM 
        staging_compound_properties
),
comparison_metrics AS (
    SELECT
        ms.source_chembl_id AS origin_molecule_id,
        ms.target_chembl_id AS comparison_molecule_id,
        mc_origin.lipophilicity_index AS origin_lipophilicity,
        mc_comparison.lipophilicity_index AS comparison_lipophilicity,
        ABS(mc_comparison.lipophilicity_index - mc_origin.lipophilicity_index) AS lipophilicity_difference
    FROM
        molecule_similarities ms
    JOIN
        molecule_characteristics mc_origin 
        ON ms.source_chembl_id = mc_origin.molecule_identifier
    JOIN
        molecule_characteristics mc_comparison 
        ON ms.target_chembl_id = mc_comparison.molecule_identifier
)
SELECT
    origin_molecule_id,
    ROUND(AVG(lipophilicity_difference), 4) AS avg_lipophilicity_difference
FROM
    comparison_metrics
GROUP BY
    origin_molecule_id;
