#' Take a clusters results and performs functional analysis enrichment for an
#' human organism
#' @param drugs
#' @param drug_prot
#' @param allowed_proteins
#' @param protein_gene
#' @param cl_id
#' @param thr
#' @param filter_type
#' @return enrichment performed
#' @authro Fernando Moreno Jabato <jabato@uma.es>
#' @import clusterProfiler enrichGO enrichKEGG
#' @import org.Hs.eg.db
#' @import KEGG.db
enrich_clusters <- function(drugs, drug_prot, allowed_proteins, 
                            protein_gene, cl_id = -1, thr = 0.001, 
                            filter_type = "None"){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  suppressPackageStartupMessages(require(KEGG.db))

  # Obtain proteins related to drug set
  proteins <- drug_prot$Protein[which(drug_prot$Drug %in% drugs)]
  proteins <- proteins[which(proteins %in% allowed_proteins)]
  # Check
  if(length(proteins) <= 0){
    return(data.frame(Cluster = numeric(0),
              Type = character(0),
              Onto = character(0),
              Pval = numeric(0),
              Term = character(0),
              Name = character(0),
              Sig_Genes = character(0)))
  }

  # Obtain genes
  genes <- protein_gene$ENTREZID[which(protein_gene$UNIPROT %in% proteins)]

  # Perform GO enrichment
  go <- clusterProfiler::enrichGO(
          gene          = genes, #genes,
          OrgDb         = "org.Hs.eg.db", #organism,
          keyType       = 'ENTREZID', #keyType,
          ont           = 'ALL', # SubOntology
          pvalueCutoff  = thr, #pvalueCutoff,
          pAdjustMethod = "BH") #qvalueCutoff)

  # Obtain info
  if(nrow(go) <= 0){
    go_df <- data.frame(Cluster = numeric(0),
              Type = character(0),
              Onto = character(0),
              Pval = numeric(0),
              Term = character(0),
              Name = character(0),
              Sig_Genes = character(0))
  }else{
    go_df <- data.frame(Cluster = rep(cl_id,nrow(go)),
              Type = rep(filter_type,nrow(go)),
              Onto = rep("GO",nrow(go)),
              Pval = go$p.adjust,
              Term = go$ID,
              Name = go$Description,
              Sig_Genes = go$geneID,
              stringsAsFactors = FALSE)
  }

  # Performs KEGG enrichment
  kegg <- clusterProfiler::enrichKEGG(
            gene         = genes, #genes,
            organism      = "hsa", #organism,
            keyType       = "kegg", #keyType,
            pvalueCutoff  = thr, #pvalueCutoff,
            pAdjustMethod = "BH", #pAdjustMethod,
            use_internal_data = TRUE)  

  kegg <- as.data.frame(kegg)

  if(nrow(kegg) <= 0){
    kegg_df <- data.frame(Cluster = numeric(0),
              Type = character(0),
              Onto = character(0),
              Pval = numeric(0),
              Term = character(0),
              Name = character(0),
              Sig_Genes = character(0))
  }else{
    kegg_df <- data.frame(Cluster = rep(cl_id,nrow(kegg)),
              Type = rep(filter_type,nrow(kegg)),
              Onto = rep("KEGG",nrow(kegg)),
              Pval = kegg$p.adjust,
              Term = kegg$ID,
              Name = kegg$Description,
              Sig_Genes = kegg$geneID,
              stringsAsFactors = FALSE)
  }

  # Collapse results
  if(nrow(go_df) <= 0 & nrow(kegg_df) <= 0){
    return(go_df)
  }else if(nrow(go_df) <= 0){
    return(kegg_df)
  }else if(nrow(kegg_df) <= 0){
    return(go_df)
  }else{
    return(rbind(go_df,kegg_df))
  }
}
