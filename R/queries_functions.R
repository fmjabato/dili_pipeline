# This function retrieves data from a relational database
# whose database-management system (DBMS) is MySQL. The information
# is retrieved executing an SQL query.
# @param username: A string containing the user name.
# @param password: A string containing the user password.
# @param host: A string containing the host direction.
# @param port: An integer containing the port number.
# @param dbname: A string containing the name of the database to be queried.
# @param query: A string containing the SQL query to be executed.
# @return a dataframe containing the results of the query, where the columns
#            are the attributes and the rows are the tuples.
# @import RMySQL dbConnect MySQL dbGetQuery dbDisconnect
# @author: Guillermo López García (guilopgar@uma.es)
# @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
mysqlQuery <- function(username, password, host, port, dbname, query) {
  con <- RMySQL::dbConnect(RMySQL::MySQL(),
                   username = username,
                   password = password,
                   host = host,
                   port = port,
                   dbname = dbname)
  
  res <- RMySQL::dbGetQuery(con, query)
  RMySQL::dbDisconnect(con)
  return(res)
}


# This function retrieves data from a SPARQL graph database
# @param endpoint: A string containing the endpoint.
# @param query: A string containing the RDF query to be executed.
# @return a dataframe containing the results of the query
# @import SPARQL SPARQL
# @author: Fernando Moreno Jabato (jabato<at>uma<dot>es)
sparqlQuery <- function(endpoint,query){
  results <- SPARQL::SPARQL(url = opt$endpoint, query = query)$results
  return(results)
}