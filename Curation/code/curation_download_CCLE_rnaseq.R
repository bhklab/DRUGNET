
## path to files
path.data <- file.path("data", "CCLE_rnaseq")
path.cell <- file.path(path.data, "celline")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=TRUE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=TRUE, recursive=TRUE) }

### Adrian please put your code to extract cell line name from the manifest here
## So the manifest of all ccle rnaseq cell lines should get download from CGhub to path.cell and then parse the html file and create a csv file named cell_line_annotation_CCLE_rnaseq.csv as output