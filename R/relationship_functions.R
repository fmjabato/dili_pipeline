# Given a Drug-Protein and a Protein-Family table, this function creates a 
# Drug-Family table, in which drugs are associated with the domains of their 
# single domain protein targets.
# The process is done by performing a database join operation on the Uniprot ID 
# column of both tables. This function is an implementation of the "PPDMs 
# mapping of small molecule bioactivities from ChEMBL to Pfam-A protein 
# domains", described in Kruger, F.A. et al. (2014).
# @param chemblTable: A data frame representing the Drug-Protein table, which 
#         must have one column containing the targets Uniprot IDs.
# @param colUniprotChembl: A positive integer representing the column of the 
#         chemblTable that contains the targets Uniprot IDs.
# @param colChembl: A positive integer or an integer vector representing the 
#         column/s of the chemblTable to be included in the returned table, 
#         e.g. the column containing the ChEMBL drugs names.
# @param pfamTable: A data frame representing the Protein-Family table, which 
#         must have one column containing the proteins Uniprot IDs.
# @param colUniprotPfam: A positive integer representing the column of the 
#         pfamTable that contains the proteins Uniprot IDs.
# @param colPfam: A positive integer or an integer vector representing the 
#         column/s of the pfamTable to be included in the returned table, e.g. 
#         the column containing the Pfam-A families names (recommended column).
# @param colSort: A positive or an integer vector representing the column/s 
#         used to optionally sort the returned Drug-Family table. If NULL, no 
#         sorting is done.
# @return a dataframe representing the Drug-Family table, containing the 
#         selected columns from the ChEMBL Drug-Protein table and then the 
#         selected ones from the Pfam-A Protein-Family table.
# @author: Guillermo López García (guilopgar@uma.es)
# @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
ppdmsDDomain <- function(chemblTable, pfamTable, colUniprotChembl = 7, 
                         colChembl = 3, colUniprotPfam = 1, colPfam = 2, 
                         colSort = c(1, 2)) {  
  # Extracting single domain proteins
  protFreqs <- table(pfamTable[[colUniprotPfam]])
  pfamSD <- pfamTable[pfamTable[[colUniprotPfam]] %in% 
                      names(protFreqs[protFreqs==1]), ]
  
  # Joining Drug-Protein and SDProtein-Family tables into a Drug-Family table
  durgDomain <- unique(merge(chemblTable[, c(colUniprotChembl, colChembl)], 
                             pfamSD[, c(colUniprotPfam, colPfam)], 
                             by.x = 1, by.y = 1)[, -1])
  
  # Sorting the table
  if(!is.null(colSort)) {
    if(!all(colSort %in% seq(length(durgDomain)))) {
      stop(paste0("The value of the 'colSort' argument must be one or more",
        " column numbers of the Drug-Family table"))
    }
    durgDomain <- durgDomain[do.call(order, lapply(durgDomain[colSort], 
                                                    as.vector)), ]
  }
  
  return(durgDomain)
}



# Given a Drug-Protein and a Protein-Family table, this function creates a 
# Drug-Family table, in which drugs are associated with the families that are 
# significantly overrepresented among their targets. The process is done by 
# performing a binomial test for every Drug-Family pair in wich the family 
# contains at least one of the targets of the drug. Then, we filter the p-values
# (optionally corrected) that correspond to an overrepresentation that is 
# statistically significant. This function is an implementation of the 
# Drug-Family mapping strategy used in Moya Garcia A., Dawson N.L., Kruger F.A.,
# et al. (2016).
# @param chemblTable: A data frame representing the Drug-Protein table, which 
#        must have one column containing the targets Uniprot IDs.
# @param colUniprotChembl: A positive integer representing the column of the 
#        chemblTable that contains the targets Uniprot IDs.
# @param colChembl: A positive integer representing the column of the 
#        chemblTable that contains the drugs identifiers to be included in the 
#        returned table, e.g. the ChEMBL names.
# @param pfamTable: A data frame representing the Protein-Family table, which 
#        must have one column containing the proteins Uniprot IDs.
# @param colUniprotPfam: A positive integer representing the column of the 
#        pfamTable that contains the proteins Uniprot IDs.
# @param colPfam: A positive integer representing the column of the pfamTable 
#        that contains the families identifiers to be included in the returned 
#        table, e.g. the Pfam-A families names.
# @param pvalue: A string representing the method used to produce the p-values 
#        when performing the binomial tests. The posible values are: "one", 
#        which means greater one-sided p-value, and "two", which means two-sided
#        p-value. If two-sided p-values are computed, the ones that represent 
#        overrepresentation are filtered. 
# @param statThres: A positive number representing the threshold used to filter 
#        the statistically significant p-values (or adjusted p-values).
# @param corrMethod: A string containing the name of the correction method used 
#        to optionally adjust the binomial tests p-values. The available methods
#        are the ones used in the p.adjust function of stats library. If NULL, 
#        no correction is done.
# @param colSort: A positive or an integer vector representing the column/s used
#        to optionally sort the returned Drug-Family table. If NULL, no sorting 
#        is done.
# @return a dataframe representing the Drug-Family table, containing 3-4 
#        columns: the drugs identifiers, the families identifiers, the binomial 
#        tests p-values and optionally the q-values (corrected p-values).
# @import stats binom.test p.adjust
# @author: Guillermo López García (guilopgar@uma.es)
# @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
overrepDFamily <- function(chemblTable, pfamTable, colUniprotChembl = 3, 
                           colChembl = 2, colUniprotPfam = 1, colPfam = 3, 
                           pvalue = "one", statThres = 0.05, corrMethod = "BH", 
                           colSort = c(1, 2)) {  
  # Join Drug-Protein and Protein-Family tables into a Drug-Family table
  drugFamily <- unique(merge(chemblTable[, c(colUniprotChembl, colChembl)], 
                                 pfamTable[, c(colUniprotPfam, colPfam)], 
                                 by.x = 1, by.y = 1))[c(2, 1, 3)]
  names(drugFamily) <- c("ChEMBL_ID", "UniProt_ID", "PfamA_ID")
  
  
  # Compute the total number of drug target proteins
  n <- length(unique(drugFamily$UniProt_ID))
  
  
  # For every protein family, compute the expected probability that any drug 
  # target is a member of that family: 
  #    Pff = the number of drug target proteins associated to a family / n
  pff <- table(unique(drugFamily[c("UniProt_ID", "PfamA_ID")])["PfamA_ID"])/n
  
  
  # For every drug, compute the number of targets (number of trials in the 
  # binomial test)
  nTarget <-table(unique(drugFamily[c("ChEMBL_ID", "UniProt_ID")])["ChEMBL_ID"])
  
  
  # For every drug-family pair, compute the number of targets of that drug that 
  # belong to that family(number of succeses in the binomial test)
  DFPvalues <- as.data.frame(table(drugFamily[c("ChEMBL_ID", "PfamA_ID")]), 
                                         stringAsFactors = FALSE)
  DFPvalues <- DFPvalues[DFPvalues$Freq > 0, ]
  
  
  # For every drug-family pair, compute either the binomial test one-sided 
  # (greater) p-value, i.e. alternative hypothesis: 
  #    p > Pff, or the two sided p-value, i.e. p != Pff
  pvalueMethod <- "g"
  if(pvalue == "two") {
    # For every drug-family pair, compute the expected number of successes 
    # (the expected number of targets of that drug that are expected to belong 
    # to that family)
    DFPvalues$e_value <- apply(DFPvalues, 1,function(x){
      nTarget[x["ChEMBL_ID"]] * pff[x["PfamA_ID"]]
    })
    
    # Compute the two sided p-value only for the overrepresented families
    DFPvalues <- DFPvalues[DFPvalues$Freq > DFPvalues$e_value, ]
    pvalueMethod <- "t"
  }
  
  DFPvalues$p_value <- apply(DFPvalues, 1, function(x) {
    round(stats::binom.test(x = as.numeric(x["Freq"]),
                            n = nTarget[x["ChEMBL_ID"]],
                            p = pff[x["PfamA_ID"]],
                            alternative = pvalueMethod)$pvalue, 
          digits = 6)
  })
  
  DFPvalues <- DFPvalues[c("ChEMBL_ID", "PfamA_ID", "p_value")]
  
  
  # For every drug-family pair, optionally adjust its p-value using the 
  # specified correction method
  if(!is.null(corrMethod)) {
    DFPvalues$adjust_p_value <- stats::p.adjust(DFPvalues$p_value, 
                                                method = corrMethod)
  }
  
  # Filter the statistically significant drug-family pairs using the last 
  # column, representing either the p-value or the adjusted ones
  colFilter <- ncol(DFPvalues)
  DFPvalues <- DFPvalues[DFPvalues[[colFilter]] < statThres, ]
  
  
  # Sort the Drug-Family table
  if(!is.null(colSort)) {
    if(!all(colSort %in% 1:length(DFPvalues))) {
      stop(paste0("The value of the 'colSort' argument must be one or more",
        " column numbers of the Drug-Family table"))
    }
    DFPvalues <- DFPvalues[do.call(order,lapply(DFPvalues[colSort],as.vector)),]
  }
  
  
  # Return the Drug-Family table
  return(DFPvalues)
}
