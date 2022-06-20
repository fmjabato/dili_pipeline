#! /usr/bin/env Rscript

# @author Fernando Moreno Jabato <jabato@uma.com>
# @description scripts to use input tables to calculate clusters

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
  optparse::make_option(c("-p", "--protein"), action="store", type="character",
              dest="dp_file", help=paste0("File with Drug-Protein network ",
                "(overrep/ppdms")),
  optparse::make_option(c("-d","--domain"), action="store",type="character",
              dest="dd_file", help="File with Drug-Domain table"),
  optparse::make_option(c("-c","--clust"), action="store_true", default = FALSE,
              dest="clusters", help="Apply clustering process"),
  optparse::make_option(c("-o", "--out"), action="store", type = "character", 
            default="./KernelAndClustering", dest="outfile",
            help=paste0("[OPTIONAL] Output file name (without extension)",
              "[Default = %default]")),
  optparse::make_option(c("-v", "--verbose"), action="store_true", 
              default=FALSE, dest="verbose",help="Activate verbose mode")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


################################################################################
## PIPELINE
################################################################################

### ### ### LOAD AND TRANSFORM ### ### ###
if(opt$verbose) message("Loading data")

# Load data
Drug_Protein <- read.table(file = opt$dp_file, sep = "\t", header = TRUE, 
                            stringsAsFactors = FALSE)
Drug_Domain  <- read.table(file = opt$dd_file, sep = "\t", header = TRUE, 
                            stringsAsFactors = FALSE)

# NOTE: SOLVE IN FUTURE !!!
dp_drugs    <- which(names(Drug_Protein) == "ChEMBL_Name")
dp_proteins <- which(names(Drug_Protein) == "UniProt_ID")
if("ChEMBL_Name" %in% names(Drug_Domain)){
  dd_drugs  <- which(names(Drug_Domain) == "ChEMBL_Name")
  dd_domain <- which(names(Drug_Domain) == "pfamA_acc") 
}else if("ChEMBL_ID" %in% names(Drug_Domain)){
  dd_drugs  <- which(names(Drug_Domain) == "ChEMBL_ID")
  dd_domain <- which(names(Drug_Domain) == "PfamA_ID")
}else{
  message("Error: imposible to find DDomains drug column")
  stop()
}

if(opt$verbose) message("Generating connectivity matrices.")

# Generate connectivity matrix from both datasets
dp_connectivity <- binary_connectivity_matrix(Drug_Protein[,c(dp_drugs,
                                                              dp_proteins)])
dd_connectivity <- binary_connectivity_matrix(Drug_Domain[,c(dd_drugs,
                                                             dd_domain)])

if(opt$verbose) message("Calculating HyI value for connectivity matrices:")

# Calculate relationship weights using HyI p-values
DD_HyI <- hyI_vectors(dd_connectivity, verbose = opt$verbose)
DP_HyI <- hyI_vectors(dp_connectivity, verbose = opt$verbose)


# Declare useful function
vector_angle <- function(x,y,module = sqrt(x^2+y^2)){
  # Check
  if(x == 0 | y == 0){
    if(module == 0) return(-1) # Both are zero
    else if(x == 0) return(90) # Vector over y-axis
    else return(0)             # Vector over x-axis
  }else{ # Both components are >0
    # Obntain alpha = angle with x-axis
    alphasin <- y/module
    alpha <- asin(alphasin)
    # Return in degrees
    return(alpha*180/pi)
  }
}

# Standarization function
standard_scale <- function(x){(x-min(x))/(max(x)-min(x))}



# Calculates networks
threshold_exps <- seq(2,8)
if(opt$verbose) message(paste0("Generating networks with several thresholds: ",
                              "10^-{",paste(threshold_exps,collapse = ","),"}"))


Edges_lists <- lapply(threshold_exps,function(x){
  # Prepare current set threshold
  threshold <- 10^(-x)
  # Apply threshold
  DD_HyI[DD_HyI > threshold] <- 1
  DP_HyI[DP_HyI > threshold] <- 1
    
  # Transform p-values to intensities and normalize
  DD_HyI <- -log(DD_HyI) # P-value to intensity
  DP_HyI <- -log(DP_HyI) 
  DD_HyI <- standard_scale(DD_HyI) # Normalize [0-1]
  DP_HyI <- standard_scale(DP_HyI) 
  
  # Obtain all drugs with any information and 
  dd_drugs <- rownames(DD_HyI)
  dp_drugs <- rownames(DP_HyI)
  drugs <- unique(c(dd_drugs,dp_drugs))
  
  # Generate all possible combinations (non_directional and without repetitions)
  combs <- expand.grid(seq_along(drugs),seq_along(drugs)) # Direc and with repe
  combs <- combs[,c(2,1)]                                 # Swap, avoid sorting
  combs <- combs[which(combs[,1] < combs[,2]),]           # Remove dir&repe
  
  # Generate a data frame of element-element values and module (relationship)
  edges <- apply(combs,1,function(indx){
    # Find drugs coords
    i_dd <- which(dd_drugs == drugs[indx[1]]) # D-d coords
    j_dd <- which(dd_drugs == drugs[indx[2]])
    i_dp <- which(dp_drugs == drugs[indx[1]]) # D-P coords
    j_dp <- which(dp_drugs == drugs[indx[2]])
    
    # Init
    dd <- 0
    dp <- 0
    
    # Take D-P info
    if(length(i_dd) != 0 & length(j_dd) != 0){ # There're info
      dd <- DD_HyI[i_dd,j_dd]
      # Check
      if(is.na(dd)) dd <- 0
    } # ELSE: There is not info
    
    # Take D-D info
    if(length(i_dp) != 0 & length(j_dp) != 0){ # There're info
      dp <- DP_HyI[i_dp,j_dp]
      if(is.na(dp)) dp <- 0
    }# ELSE: There is not info
    
    # Check
    if(dd == 0 & dp == 0){
      return(NULL)
    }
    
    # Obtain module value
    module <- sqrt(dp^2+dd^2)
    
    # Angle
    angle <- vector_angle(dp,dd)
    
    # Store and return
    return(list(nodeA  = drugs[indx[1]],
                nodeB  = drugs[indx[2]],
                dd     = dd,
                dp     = dp,
                module = module,
                angle  = angle))
  })
  
  # Remove NULLs
  edges <- edges[-which(unlist(lapply(edges,function(item){
    return(is.null(item))
  })))]
  
  # Transform edges in a dataframe
  edges <- as.data.frame(do.call(rbind,edges),stringsAsFactors = F)
  invisible(lapply(seq_along(names(edges)), function(i){
    edges[,i] <<- unlist(edges[,i])
  }))
  
  return(edges)
})




# Obtain number of nodes per cut
num_nodes <- c(length(unique(c(Drug_Protein[,dp_drugs], 
                               Drug_Domain[,dd_drugs]))),
               unlist(lapply(Edges_lists,function(edges){
                return(length(unique(c(edges$nodeA,edges$nodeB))))
              })))
names(num_nodes) <- c("original",paste("10^-",threshold_exps,sep=""))

### ### ### FIRST STORE ### ### ###
############ Verbose point
outF <- paste(opt$outfile,".RData")
if(opt$verbose) message(paste0("Storing Edges list into: ", outF))

# Store
save(Edges_lists,file = outF)


### ### ### CLUSTERING ### ### ###
if(opt$clusters){
  if(opt$verbose) message(paste0("Clustering network (threshold 10^-3). It can",
      " take a while:"))
  
  # Obtaing clusters
  LinksTime <- system.time(Links <- linkcomm::getLinkCommunities(
    Edges_lists[[2]][,c("nodeA","nodeB","module")],
    plot = F, verbose = opt$verbose))
  
  if(opt$verbose){
    message(paste("Clustering finished. Time:",LinksTime[["elapsed"]],"s"))
    message("Storing clustering with edges lists already stored")
  }
  # remove unncessary variables
  attach(list("outfile"=opt$outfile,"verbose"=opt$verbose))
  wanted <- c("num_nodes","Edges_lists","Links","dp_connectivity","dd_connectivity","opt")
  rm(list = ls()[-which(ls() %in% wanted)])
  # Store
  save.image(file = paste0(outfile,".RData"))
}

