# Function to obtain PFam domain name using Pfam REST API.
# Usage of this functiuon is internal. Inputs will not be checked
# @param id: of domain family wanted
# @param endpoint: of REST API
# @return name of famyly asked
# @import RCurl getURL
# @authro Fernando Moreno Jabato <jabato@uma.es>
pfam_domain_family_name <- function(id, 
                                    endpoint = "http://pfam.xfam.org/family/"){  
  xml_doc <- RCurl::getURL(url = paste0("http://pfam.xfam.org/family/",
                                        id,"?output=xml"))
  name <- gsub(".*<description>\n<!\\[CDATA\\[\n","",xml_doc)
  name <- gsub("\n.*","",name)
  return(name)
}
