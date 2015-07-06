library(shiny)

shinyUI(
  
  fluidPage(
      titlePanel("Drug Similarity Search"),
      sidebarLayout(
          sidebarPanel(
              selectInput("optn", label = "Search by:", choice = c("SMILES String", "Chemical Name")),
              uiOutput("prompt"),
              checkboxGroupInput("studies", label = h4("Studies"), list("CCLE", "CGP", "GRAY", "GNE", "GSK", "CTRP", "NCI60")),
              sliderInput("cutoff", label = h4("Cutoff in %"), min = 0, max = 100, value = 50),
              checkboxInput("row", label = "Display one row per drug", value = FALSE),
              checkboxInput("smiles", label = "Display SMILES", value = FALSE),
              checkboxInput("inchikey", label = "Display InChIKey", value = FALSE),
              checkboxInput("chembl", label = "Display links to ChEMBL", value = FALSE),
              actionButton("submitButton", "Submit!")
              ),
          mainPanel(
              p("Search for similar drugs in large-scale pharmacogenomic studies to a given drug according using the Tanimoto algorithm."),
              p("Current studies include CCLE, CGP, GRAY, GNE, GSK, CTRP, and NCI60."),
              p("By default, an input drug will be searched against drugs in all studies."),
              textOutput("params"),
              h3("Results"),
              textOutput("record"),
              textOutput("res"),
              dataTableOutput("similarDrugs")
        )
      )
      )
  )

