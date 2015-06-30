library(shiny)
library(rcdk)
library(webchem)
library(RJSONIO)
library(RCurl)

## Global variables
load("data_obj.RData")
load("fingerprints.RData")



shinyServer(function(input, output) {
    
  output$prompt <- renderUI({
        switch(input$optn,
        "smiles" = textInput("string", label = h4("Input SMILES string")),
        "name" = textInput("string", label = h4("Input chemical name"))
        )
    })
  
  output$params <- renderText({
      paste("You searched for all drugs matching", input$string, "with at least", input$cutoff, "% similarity in",
            paste0(input$studies, collapse = ","))
    })

  output$similarDrugs <- renderDataTable({

    toConvert <- switch(input$optn, "smiles" = input$string, "name" = getSMILESfromName(input$string))
    query.mol <- tryCatch(parse.smiles(toConvert)[[1]], error = function (e) { renderText({"No SMILES string entered"}) })
    vals <- findSimilarDrugs(query.mol, input$cutoff, input$studies)
        
    output$res <- renderText({
       paste(nrow(vals), "records returned.")
     })
    
  # if(!input$display) {
     # vals <-cbind(vals, appendIdentifers(vals, input$studies))
    #}
    
   
   # if (input$chembl) {
     # names <- as.vector(t(vals["names"]))
      #chemblIDs <- lapply(names, getChEMBLfromName)
      #links <- createCHEMBLPages(chemblIDs, names)
      #vals["names"] <- links
   #}
  #  rownames(vals) <- NULL
    return(vals)
  }, escape = FALSE)

})

