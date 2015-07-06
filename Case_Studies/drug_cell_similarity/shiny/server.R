library(shiny)
library(rcdk)

source("chemicalTransUtils.R") ## Functions to translate between chemical identifers
load("drug_dat_with_fps.RData") ## database with drug curation & corresponding ids/fingerprints
load("chembls.RData") ## Load chembl Ids

## Use tanimoto algorithm to return a data.frame of similar drugs against all drugs in DB
## smiles: a smiles id 
findSimilarDrugs <- function(smiles) {
  
  ## Get fingerprint for drug and compare it amongst all drugs in search space
  query.fp <- get.fingerprint(smiles, type="extended")
  sims <- unlist(lapply(drugs$fps, distance, fp2 = query.fp, method = "tanimoto"))
  
  ## Get unique entries from frame and order
  frame <- data.frame("drug_id" = drugs$drug_id, "smiles" = drugs$smiles, 
                      "inchikey"= drugs$inchikey, "study_name" = drugs$study_name, 
                      "similarity_score" = sims)
  frame <- frame[order(frame$similarity_score, decreasing = TRUE), ]
  
  return(frame)
}

## Create hyperlink to chembl record around chemical name, if a chembl record exists
## name: a chemical name
createCHEMBLPage <- function(name) {
    if (is.null(name)) {
      "No ChEMBL ID Found" 
    }
    else { ## Create link for each unique chembl ID
      paste0(
      lapply(name, function (x) sprintf('<a href="https://www.ebi.ac.uk/chembl/compound/inspect/%s">%s</a>', x, x)),
      collapse = ","
    )
    }
}

## Create hyperlinks based on all names
## name: vector of chemical name
createCHEMBLPages <- function(drugs) {
  links <- sapply(chembls[as.character(drugs)], createCHEMBLPage) ## Get ChEMBL IDs from chembl list then create pages
  return(links)
}

shinyServer(function(input, output) {
  
  output$prompt <- renderUI({
    textInput("string", label = h4(paste("Input", input$optn))) 
  })
    
  output$params <- renderText({
    
    if (length(input$studies) == 0) {
      paste("You searched for all drugs matching", input$string, "with at least", input$cutoff, "% similarity in all
             available studies.")
    }
    
    else {
      paste("You searched for all drugs matching", input$string, "with at least", input$cutoff, "% similarity in",
            paste0(input$studies, collapse = ","))
    }
    
    })

  output$similarDrugs <- renderDataTable({ 
    
    input$submitButton
    
    isolate({
    toConvert <- switch(input$optn, "SMILES String" = input$string, "Chemical Name" = getSMILESfromName(input$string))
    query.mol <- parse.smiles(toConvert)[[1]]
    vals <- findSimilarDrugs(query.mol)
    })
    
    vals <- vals[vals$similarity_score > (input$cutoff / 100), ] ## Cutoff according to use defined score
    
    if (length(input$studies) != 0) {
      vals <- vals[which(vals$study_name %in% input$studies), ] ## Cutoff according to studies
    }
    
    if (input$row) {
      ## Return one row for drug (i.e. remove duplicates)
      vals$studies <- sapply(vals$drug_id, 
                              function (name) paste0(t(vals[which(vals$drug_id == name), "study_name"]), collapse = ","))
      vals <- vals[!duplicated(vals$drug_id), -which(names(vals) %in% "study_name")]
    }
    
    output$record <- renderText({
      paste(nrow(vals), "records returned.")
    })
    
    if(!input$inchikey) {
     vals <- vals[, -which(names(vals) %in% c("inchikey"))]
    }
    
    if (!input$smiles) {
      vals <- vals[, -which(names(vals) %in% c("smiles"))]
    }
    
    if (input$chembl) {
       vals$chembls <- createCHEMBLPages(vals$drug_id)
    }
    
    return(vals)
  }, escape = FALSE) 

})

