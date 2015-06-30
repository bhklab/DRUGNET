## Load files
load("Expression.arrays.Rdata")
load("merged_annotations.RData")

## Process SNV data
## SU2C is a subset of gray data.
snps <- data.frame("unique.cellid" = cell.line.anno$unique.cellid, "ccle" = cell.line.anno$ccle.filename, "cgp" = cell.line.anno$cgp.filename, "gray" = cell.line.anno$su2c.filename)

ccleExpDat <- data.frame("unique.cellid" = sampleinfo.ccle$cellid, "name" = sampleinfo.ccle$filename, "array" = sampleinfo.ccle$Expression.arrays)

cgpExpDat <- data.frame("unique.cellid" = sampleinfo.cgp$cellid, "name" = sampleinfo.cgp$filename, "array" =  sampleinfo.cgp$Array.Data.File)


snpStudies <- c("CCLE", "CGP", "GRAY")
types <- c("gene_expression_array", "SNP")

importGeneExpressionType <- function(type) {
	connection <- connDB()
	dbGetQuery(connection, paste0("INSERT INTO GENE_EXPRESSION_TYPE (stat_name) VALUES (\'", type, "\')"))
	dbDisconnect(connection)
}

getGeneExpressionType <- function(type) {
	connection <- connDB()
	id <- dbGetQuery(connection, paste0("SELECT stat_id FROM GENE_EXPRESSION_TYPE WHERE stat_name = \'", type, "\'"))
	dbDisconnect(connection)
	return(id)
}
sapply(types, importGeneExpressionType)

writeSNP <- function(study) {
	frame <- data.frame(
	"cellline_id" = snps$unique.cellid,
	"exp_study_id" = getStudyID(study),
	"stat_id" = getGeneExpressionType('snp'),
	"file_path" = snps[, tolower(study)]
	)
	frame <- frame[complete.cases(frame),]
	con <- connDB()
	dbWriteTable(con, "GENE_EXPRESSION_METADATA", frame, append = TRUE, row.names = FALSE)
	dbDisconnect(con)
}

sapply(snpStudies, writeSNP)

writeGene <- function(study, obj) {
	frame <- data.frame(
		"cellline_id" = obj$unique.cellid,
		"exp_study_id" = getStudyID(study),
		"stat_id" = getGeneExpressionType("gene_expression_array"),
		"file_path" = obj$name)
	frame <- frame[complete.cases(frame),]	
	con <- connDB()
	dbWriteTable(con, "GENE_EXPRESSION_METADATA", frame, append = TRUE, row.names = FALSE)
	dbDisconnect(con)
}

mapply(writeGene, c("CCLE", "CGP"), list(ccleExpDat, cgpExpDat))
