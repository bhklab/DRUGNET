### query = "provide a list of queries in the form of \"field:value\" to match your results against." ###
### fields = "List of fields will determine the fields returned" ###
### sort = "Field by which to sort" ###
### descending = "Sort descending or ascending" ####
### skip = "How many lines to skip" ###
### limit = "How many entries to return" ###
### returnMatrix = "Return a matrix of the pulled values, with multiple entries concatecated by collapse. If false, returns a nested list representing the JSON object."

#user_key="afe4848514df2b491be1c49b49475cab"



lincsAPIQuery <- function(infoType = c("pertinfo", "cellinfo", "geneinfo", "instinfo", "plateinfo", "siginfo"), d=NULL, query=NULL, fields=NULL, sorted=NULL, descending = FALSE, skip=0, limit=999999999, returnMatrix=TRUE, user_key=stop("Please Provide a LINCS api key")) {
  
  require(RCurl)
  require(rjson)
  
  infoType <- match.arg(infoType)
  
  urlString <- paste("/a2/", infoType, "?", sep="")
  
  queryString <- "q={"
  
  
  
  if (!is.null(query)){
    
    pieces <- unlist(strsplit(query[1], split = ":"))
    queryString <- paste(queryString, "\"", pieces[1], "\":\"", pieces[2], "\"", sep="")
    
    for (q in query[-1]){
      pieces <- unlist(strsplit(q, split = ":"))
      queryString <- paste(queryString, ",\"", pieces[1], "\":\"", pieces[2], "\"", sep="")
      
    }
  }
  
  queryString <- paste(queryString, "}", sep="")
  
  urlString <- paste(urlString, queryString, sep="")
  
  if (!is.null(fields)){
    
    fieldsString <- "&f={"
    
    fieldsString <- paste(fieldsString, "\"", fields[1], "\":1\"", sep="")
    
    for (f in fields[-1]){
      
      fieldsString <- paste(fieldsString, "\"", f[1], "\":1\"", sep="")
      
    }
    fieldsString <- paste(fieldsString, "}", sep="")
    urlString <- paste(urlString, fieldsString, sep="")
  }
  
  if (!is.null(sorted)){
    
    ###### 1 for ascending; -1 for descending
    
    if(descending){
      
      sortString <- paste("&s={\"",sorted,"\":-1}", sep="")
      
    } else {
      sortString <- paste("&s={\"",sorted,"\":1}", sep="" )
    }
    
    urlString <- paste(urlString, sortString, sep="")
    
  }
  
  if (!is.null(d)){
    
    dString <- paste("&d=", d, sep="")
    
    urlString <- paste(urlString, dString, sep="")
    returnMatrix=FALSE
  }
  
  if (skip>0) {
    
    skipString <- paste("&sk=", skip, sep="")
    urlString <- paste(urlString, skipString, sep="")
  }
  
  urlString <- paste(urlString, "&l=", limit, sep="")
  urlString <- paste(urlString, "&user_key=", user_key, sep="")
  
  
  response <- getURI(url = paste("http://api.lincscloud.org", urlString, sep=""))
  
  result <- fromJSON(json_str = response)
  if (!returnMatrix){ return(result)} else {
    dm <- matrix(nrow=length(result), ncol=length(unique(unlist(sapply(result, names)))))
    
    colnames(dm) <- unique(unlist(sapply(result, names)))
    
    for (i in 1:length(result)){
      for (name in names(result[[i]])){
        dm[i,name] <- paste(result[[i]][[name]], collapse="///")
      }
    }
    
    return(dm)
  }
}