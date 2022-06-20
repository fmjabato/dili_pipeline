#! /usr/bin/env Rscript

#' @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
#' @description: Combines Chembl and pfam tables to obtains useful triplets table

################################################################################
## CONFIGURE SCRIPT
################################################################################
# Prepare input commands
option_list <- list(
  optparse::make_option(c("-c", "--chembl"), action="store", type="character",
			  dest="chemblF", help="ChEMBL table"),
  optparse::make_option(c("-p", "--pfam"), action="store", type="character",
			  dest="pfamF", help="PFAM table"),
  optparse::make_option(c("-o", "--out"), action="store", type = "character",
			  dest="outF", help="Output file")
)

# Parse command line options
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


################################################################################
## PIPELINE
################################################################################

# Load
DP <- read.table(file = opt$chemblF, sep = "\t", quote = "", 
				stringsAsFactors = FALSE, header = TRUE)
Pd <- read.table(file = opt$pfamF, sep = "\t", quote = "", 
				stringsAsFactors = FALSE, header = TRUE)



triplets <- as.data.frame(do.call(rbind, lapply(unique(DP$ChEMBL_Name), 
	function(drug){
  # Find proteins into each drug
  prots <- DP$UniProt_ID[which(DP$ChEMBL_Name == drug)]
  # per each prot, obtain it domain
  info <- as.data.frame(do.call(rbind,lapply(prots,function(prot){
    # Find domains
    indx <- which(Pd$uniprot_acc == prot)
    # Check
    if(length(indx) == 0){
      return(data.frame(Drug = drug, 
                        Protein = prot, 
                        Domain_ID = NA, 
                        Domain_Name = NA, 
                        stringsAsFactors = FALSE))
    }
    # Prepare data
    prot_domains <- data.frame(Drug = rep(drug,length(indx)),
                               Protein     = rep(prot, length(indx)),
                               Domain_ID   = Pd$pfamA_acc[indx],
                               Domain_Name = Pd$pfamA_id[indx], 
                               stringsAsFactors = FALSE)
		return(prot_domains)
	})))
  return(info)
})))

# Write as csv
write.csv2(triplets,file = opt$outF, quote = FALSE, row.names = FALSE)

