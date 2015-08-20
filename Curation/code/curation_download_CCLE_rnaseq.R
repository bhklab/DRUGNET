<<<<<<< HEAD
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
=======

## path to files
path.data <- file.path("data", "CCLE_rnaseq")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=TRUE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=TRUE, recursive=TRUE) }

### Adrian please put your code to extract cell line name from the manifest here
## So the manifest of all ccle rnaseq cell lines should get download from CGhub to path.cell and then parse the html file and create a csv file named cell_line_annotation_CCLE_rnaseq.csv as output
>>>>>>> 418e467c61132c0415d1d941145733a31cae7061
