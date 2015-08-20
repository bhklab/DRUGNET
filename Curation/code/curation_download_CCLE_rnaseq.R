#### R Script to parse CCLE RNA-Seq Manifest File
#### Manifest file cannot be downloaded automatically. It must be downloaded
#### manually from https://browser.cghub.ucsc.edu/search/?study=("Other_Sequencing_Multiisolate)
#### &state=(live)&limit=15&library_strategy=(RNA-Seq) by adding all of the items
#### to cart and then downloading the manifest from that.

library(XML)

## x: name of XML file
parseDoc <- function(x) {
	
	## Obtain document tree
	doc <- xmlTreeParse(x)
	
	## Obtain root of xml tree
	top <- xmlRoot(doc)
		
	## Parse results for each result
	frames <- lapply(2:(length(top) - 1), function (x) {
		val <- xmlToList(top[[x]])
		files <- unlist(val$files) ## files needs extra parsing because it is a sub XML node
		impInfo <- data.frame(
		"analysis_id" = val$analysis_id,
		"file_name" = files[1],
		"size" = files[2],
		"analysis_data_uri"= val$analysis_data_uri) ## construct frame from important information
	})
	res <- do.call("rbind", frames)
	return(res)
	
}

## Assume in dir where all manifest files present
if (!file.exists("manifest.xml")) {
	stop("The manifest for CCLE RNA-Seq does not exist! Please download it from CGHub.")
}

## parse file into a df
frame <- parseDoc("manifest.xml")

## Create column cells with cell line name for each file
cells <- sapply(as.vector(t(frame["file_name"])), function (x) {
		str <- sub("\\.[0-9]\\.bam", "", x) ## remove end of file name
		str <- sub("\\w*\\.", "", str) ## remove front of filename
	 })

## Construct and write final frame
frame <- data.frame(frame, "cells" = cells, "tissueid" = NA)
write.csv(frame, file = "ccle_rnaseq_curation.csv")