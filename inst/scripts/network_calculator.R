#! /usr/bin/env Rscript

# @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
# @author: Guillermo López García (guilopgar@uma.es)
# @description: Drugs weighted networks calculator


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

# Prepare input commands
option_list <- list(
  optparse::make_option(c("-c", "--chembl-table"), action="store", 
              type="character", dest="chembl_file", 
              help=paste0("Absolute path of the file containing the ChEMBL",
              " Drug-Protein table in tab format with a header line")),
  optparse::make_option(c("-u", "--chembl-uniprot"), action="store", 
              type="integer", help=paste0("Column number of the ChEMBL",
                " Drug-Protein table that contains the protein targets ",
                "Uniprot IDs")),
  optparse::make_option(c("-i", "--chembl-id"), action="store", 
              type = "integer", help=paste0("Column number of the ChEMBL ",
                "Drug-Protein table that contains the drugs identifiers to be",
                " included in the Drug-Family table")),
  optparse::make_option(c("-p", "--pfam-table"), action="store", 
              type="character", dest="pfam_file", 
              help=paste0("Absolute path of the file containing the Pfam-A ",
                "Protein-Family table in tab format with a header line")),
  optparse::make_option(c("-n", "--pfam-uniprot"), action="store", 
              type="integer", help=paste0("Column number of the Pfam-A ",
                "Protein-Family table that contains the proteins Uniprot IDs")),
  optparse::make_option(c("-f", "--pfam-id"), action="store", type = "integer",
              help=paste0("Column number of the Pfam-A Protein-Family table ",
                "that contains the families identifiers to be included in the",
                " Drug-Family table")),
  optparse::make_option(c("-O", "--overrep"), action="store_true", 
              default=FALSE, help="Activate Overrepresentation network method"),
  optparse::make_option(c("-P", "--PPDMS"), action="store_true", 
              default=FALSE, help="Activate PPDMS network method"),
  optparse::make_option(c("-t", "--stat-thres"), action="store", type ="double",
              dest="stat_thres",help=paste0("[ONLY OVERREP] Threshold used to",
                " filter the statistically significant p-values (or q-values,",
                " if p-values correction is performed)")),
  optparse::make_option(c("-m", "--corr-method"), action="store", 
              type = "character", dest="corr_method",
              help=paste0("[ONLY OVERREP][OPTIONAL] Correction method used to",
              " adjust the binomial tests p-values")),
  optparse::make_option(c("-s", "--sort"), action="store", type = "integer",
              dest="column_sort",
              help=paste(paste0("[OPTIONAL] An integer (between 1 and 4) ",
              "representing the column/s used to sort the Drug-Family table:"),
        "\t\t- 1 means that the first column (ChEMBL ID) is used",
        "\t\t- 2 means that the second column (Pfam-A family ID) is used",
        "\t\t- 3 means that the third column (p-value) is used", 
        paste0("\t\t- 4 means that the fourth column (q-value, if p-values ",
          "correction is performed) is used"),
        paste0("\t\t- 5 means that the first and the second columns are used",
          " (in that particular order)"),
        paste0("\t\t- 6 means that the second and the first columns are used",
          " (in that particular order)"),
              sep = "\n")),
  optparse::make_option(c("-o", "--out"), action="store", type = "character",
              dest="out_file", help=paste0("Absolute path of the output file",
                " to be created containing the network table")),
  optparse::make_option(c("-v", "--verbose"), action="store_true", 
              default=FALSE, help="[OPTIONAL] Activate verbose mode")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


################################################################################
## PIPELINE
################################################################################

# Reading the ChEMBL Drug-Protein table
if(opt$verbose) message(paste("Loading ChEMBL Drug-Protein table from:", 
                              opt$chembl_file, "..."))
chembl <- read.table(file = opt$chembl_file, 
                     sep = "\t", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# Reading the Pfam-A Protein-Family table
if(opt$verbose) message(paste("Loading Pfam-A Protein-Family table from:", 
                              opt$pfam_file, "..."))
pfam <- read.table(file = opt$pfam_file, 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

# Sorting trick
if(!is.null(opt$column_sort) && opt$column_sort > 4) {
  aux <- c(1,2)
  opt$column_sort <- ifelse(opt$column_sort == 5,aux,rev(aux))
}


if(opt$verbose) message("Mapping ChEMBL drugs to Pfam-A protein families ...")

######################
# OVERREPRESENTATION METHOD
if(opt$overrep){
  if(opt$verbose) message("Launching Overrepresentation method")
  overrepT <- overrepDFamily(chemblTable = chembl, 
                             colUniprotChembl = opt$chembl_uniprot, 
                             colChembl = opt$chembl_id, 
                             pfamTable = pfam, 
                             colUniprotPfam = opt$pfam_uniprot, 
                             colPfam = opt$pfam_id, 
                             statThres = opt$stat_thres, 
                             corrMethod = opt$corr_method, 
                             colSort = opt$column_sort)
  if(!nrow(overrepT)) {
    message("Warning: Overrep results are empty. NO OUTPUT will be written")
  }else{
    overrepOutF <- paste0(opt$out_file,"_overrep.tsv")
    write.table(overrepT, file = overrepOutF, 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}



######################
# PPDMS METHOD
if(opt$PPDMS){
  if(opt$verbose) message("Launching PPDMS method")
  ppdmsT <- ppdmsDDomain(chemblTable = chembl, 
                         colUniprotChembl = opt$chembl_uniprot, 
                         colChembl = opt$chembl_id, 
                         pfamTable = pfam, 
                         colUniprotPfam = opt$pfam_uniprot, 
                         colPfam = opt$pfam_id, 
                         colSort = opt$column_sort)

  if(!nrow(ppdmsT)) {
    message("Warning: Overrep results are empty. NO OUTPUT will be written")
  }else{
    ppdmsOutF <- paste0(opt$out_file,"_ppdms.tsv")
    write.table(ppdmsT, file = ppdmsOutF, 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}






