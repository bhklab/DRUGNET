## Global variables 
library(RCurl)
library(RJSONIO)

## host of the web service
host <- "http://142.1.33.22:5000" 

## json:  json string containing  data
## return data.frame with parsed json string 
processJSON <- function(json, attr) {
  df <- fromJSON(json, nullValue = NA)
  df <- df[[attr]]
  return(df)
}

