/* Query syntax for a local Pfam MySQL database. This query retrieves
 * proteins and their Pfam A domains
 *
 *
 * Edited Aurelio Moya Garcia 05/07/2017 from David Velasco code
 * Edited Guillermo Lopez Garcia 10/04/18
 */

select distinct up.uniprot_acc, pf.pfamA_acc, pf.pfamA_id
FROM uniprot_reg_full up
JOIN pfamA pf ON pf.pfamA_acc = up.pfamA_acc
JOIN uniprot upr ON up.uniprot_acc = upr.uniprot_acc
WHERE upr.ncbi_taxid = '9606'
ORDER BY up.uniprot_acc, pf.pfamA_acc
