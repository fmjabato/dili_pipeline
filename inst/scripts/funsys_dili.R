#! /usr/bin/env Rscript

#' @description calculates clusters functional systems enrichments. This script
#'   has not being coded to be used in general problems. It's designed to offer
#'   an ad oc solution 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>

################################################################################
## CONFIGURE SCRIPT
################################################################################
# Obtain this script directory
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
              error=function(e) # works when using R CMD
                normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                commandArgs())], '='))[2]))
main_path_script <- dirname(full.fpath) # <root>/inst/scripts
root_path <- file.path(main_path_script, '..', '..')
curr_dir <- getwd()
libs_dir <- file.path(root_path,"R")
# Load custom libraries
custom_libraries <- list.files(path = libs_dir, pattern = "R$")
for (lib in custom_libraries){
  source(file.path(libs_dir, lib))
}


suppressPackageStartupMessages(require(org.Hs.eg.db))


# Prepare input commands
option_list <- list(
  optparse::make_option(c("-s", "--session"), action="store", type="character",
              dest="session", help="DILI results session"),
  optparse::make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file basename")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


##############################################################################
##                             LOAD & TRANSFORM                             ##
##############################################################################

# Load R Data session
load(opt$session)

# Obtain all Drug-Protein tuples
drug_prot <- as.data.frame(do.call(rbind,lapply(seq(nrow(proteins_metainfo)),
  function(i){
  # Obtain drugs
  drugs <- unlist(strsplit(proteins_metainfo[i,8],";"))
  # Return info
  return(data.frame(Protein = rep(proteins_metainfo$Protein[i],length(drugs)), 
                    Drug = drugs, stringsAsFactors = FALSE))
})))

all_drugs <- unique(drug_prot$Drug)

# Obtain target proteins
full_dili_proteins <- proteins_metainfo$Protein[
                                    which(proteins_metainfo$DILI_Score == 1)]
relevant_dili_proteins <- proteins_metainfo$Protein[
                                  which(proteins_metainfo$DILI_Score >= 0.75)]
all_proteins <- unique(proteins_metainfo$Protein)

# Obtain clusters info
clusters_nodes <- Links$nodeclusters
# Concat Cluster 0 
cl0 <- data.frame(node = all_drugs, cluster = rep(0,length(all_drugs)))
clusters_nodes <- rbind(cl0,clusters_nodes)

# Obtain Protein-Gene links
uniprots <- Rkeys(org.Hs.egUNIPROT)
protein_gene <- select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")


##############################################################################
##                                   ENRICH                                 ##
##############################################################################


# Per each cluster, obtain enrichments
clenrichments <- as.data.frame(do.call(rbind,
  lapply(unique(clusters_nodes$cluster), function(cl){
  # Enrich using only DILI proteins
  enr_dili_full <- enrich_clusters(
            drugs = clusters_nodes$node[which(clusters_nodes$cluster == cl)],
            drug_prot = drug_prot,
            allowed_proteins = full_dili_proteins,
            protein_gene = protein_gene,
            cl_id = cl, thr = 0.001, filter_type = "DILI_Full")
  # Enrich with significant DILI proteins
  enr_dili_sig <- enrich_clusters(
            drugs = clusters_nodes$node[which(clusters_nodes$cluster == cl)],
            drug_prot = drug_prot,
            allowed_proteins = relevant_dili_proteins,
            protein_gene = protein_gene,
            cl_id = cl, thr = 0.001, filter_type = "DILI_Sig")
  # Enrich with all proteins
  enr_all <- enrich_clusters(
            drugs = clusters_nodes$node[which(clusters_nodes$cluster == cl)],
            drug_prot = drug_prot,
            allowed_proteins = all_proteins,
            protein_gene = protein_gene,
            cl_id = cl, thr = 0.001, filter_type = "All")

  # Concat
  final_df <- rbind(enr_dili_full,enr_dili_sig)
  final_df <- rbind(final_df,enr_all)

  return(final_df)
})))

# Transform genes info into proteins info
clenrichments$Sig_Proteins <- unlist(lapply(
  seq_along(clenrichments$Sig_Genes), function(i){
  # Obtain cluster drugs
  cl_drugs <- clusters_nodes$node[
                    which(clusters_nodes$cluster == clenrichments$Cluster[i])]
  # Obtain related proteins
  cl_prots <- drug_prot$Protein[which(drug_prot$Drug %in% cl_drugs)]
  # Split
  genes <- unlist(strsplit(clenrichments$Sig_Genes[i],"/"))
  # Translate
  proteins <- unlist(lapply(genes,function(gene){
    # Find
    protein <- protein_gene$UNIPROT[which(protein_gene$ENTREZID == gene)]
    return(protein)
  }))
  proteins <- unique(proteins)
  proteins <- proteins[which(proteins %in% cl_prots)]
  # Return
  return(paste(proteins,collapse = "/"))
}))

##############################################################################
##                                   STORE                                  ##
##############################################################################
write.table(clenrichments, file = paste(opt$output,"_enr.tab",sep=""), 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(drug_prot, file = paste(opt$output,"_DP.tab",sep=""),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(clusters_nodes, file = paste(opt$output,"_cls.tab",sep=""),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(drugs, file = paste(opt$output,"_drugs.tab",sep=""),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(proteins_metainfo, file = paste(opt$output,"_prots.tab",sep=""),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)