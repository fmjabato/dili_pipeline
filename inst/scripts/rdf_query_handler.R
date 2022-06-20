#! /usr/bin/env Rscript

#' @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
#' @author: Guillermo López García (guilopgar@uma.es)
#' @description: RDF queries handler.


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
  optparse::make_option(c("-q", "--query"), action="store", type="character",
              dest="query_file", help=paste0("Absolute path of the file ",
                "containing the SPARQL query to be executed")),
  optparse::make_option(c("-e", "--endpoint"), action="store", type="character",
              help="URI of the SPARQL endpoint service to be queried"),
  optparse::make_option(c("-t", "--targets"), action="store", type ="character",
              help=paste0("[OPTIONAL] Absolute path of the text file ",
                "containing the target IDs to be filtered, separated either",
                " by a space or by an end of line")),
  optparse::make_option(c("-c", "--column"), action="store", type = "integer",
              help=paste0("[OPTIONAL] The number of the column of the query",
                " results that contain the target IDs to be filtered")),
  optparse::make_option(c("-o", "--out"), action="store", type = "character",
              dest="out_file", help=paste0("Absolute path of the output file ",
                "to be created with the results of the query")),
  optparse::make_option(c("-v", "--verbose"), action="store_true",default=FALSE,
              help="Activate verbose mode")
)

# Parse command line options
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


################################################################################
## PIPELINE
################################################################################

# Executing the query
if(opt$verbose) message("Executing query")
query <- readChar(opt$query_file, nchars = file.info(opt$query_file)$size)
queryResult <- sparqlQuery(endpoint = opt$endpoint, query = query)

# Check whether the query has produced results
if(!nrow(queryResult)) {
  message("Warning: The query answer is empty. NO OUTPUT will be written")
}else{
  # Filtering the target IDs
  if(!is.null(opt$targets)) {
    if(opt$verbose) message(paste("Filtering IDs by:", opt$targets))
    targetID <- scan(opt$targets, 
                      what = typeof(queryResult[[opt$column]]), 
                      sep = "\n", quiet = TRUE)
    targetsIndx <- which(toupper(queryResult[[opt$column]]) %in% 
                         toupper(targetID))
    if(length(targetsIndx)>0) queryResult <- queryResult[targetsIndx,]
  }

  if(opt$verbose) message(paste("Exporting results to file:", opt$out_file))
  write.table(queryResult, 
              file = opt$out_file, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
}

if(opt$verbose) message("Query pipeline finished.")