/* Query syntax for a local ChEMBL database. This query retrieves
 * a subset of 'druggable' targets and their ligands. It uses JOIN to
 * reduce the time complexity.
 * 
 * Search criteria is described in Moya Garcia A., Dawson N.L., Kruger F.A., et al. (2016). http://doi.org/10.1038/s41598-017-10012-x
 *
 * Author: Aurelio Moya Garcia <amoyag@uma.es>
 *
 * Edited by: Guillermo Lopez Garcia <guilopgar@uma.es>
*/

SELECT DISTINCT md.CHEMBL_ID AS 'ChEMBL_ID', 
                md.pref_name AS 'ChEMBL_Name', 
                cs.accession AS 'UniProt_ID'

FROM molecule_dictionary                md
        JOIN activities                 act USING (molregno)
        JOIN assays                     ass USING (assay_id)
        JOIN target_dictionary          td USING (tid)
        JOIN target_components          tc USING (tid)
        JOIN component_sequences        cs USING (component_id)

WHERE act.pchembl_value >= 5
        AND (data_validity_comment IS NULL OR data_validity_comment = 'manually validated')
        AND md.max_phase = 4
        AND assay_type IN ('B','F')
        AND relationship_type = 'D'
        AND target_type = 'SINGLE PROTEIN'
        AND standard_relation = '='
        AND cs.organism = "Homo sapiens"

ORDER BY md.CHEMBL_ID, md.pref_name, cs.accession