library(shiny)

shinyUI(fluidPage(
  titlePanel("Drug Similarity Search"),
  sidebarLayout(
    sidebarPanel(
      p("Search by:"),
      radioButtons("optn", label = "", choice = list("SMILES String" = "smiles", "Chemical Name" = "name")),
      uiOutput("prompt"),
      checkboxGroupInput("studies", label = h4("Studies"), list("CCLE", "CGP", "GRAY", "GNE", "GSK", "CTRP", "NCI60", "all studies")),
      sliderInput("cutoff", label = h4("Cutoff"), min = 0, max = 100, value = 50),
      checkboxInput("display", label = "Display SMILES and InChIKeys of Drugs", value = FALSE),
      checkboxInput("chembl", label = "Display links to ChEMBL", value = FALSE),
      submitButton("Submit!")),
    mainPanel(
      p("Search for similar drugs in large-scale pharmacogenomic studies to a given drug according using the Tanimoto algorithm."),
      p("Current studies include CCLE, CGP, GRAY, GNE, GSK, CTRP, and NCI60"),
      textOutput("params"),
      h3("Results"),
      textOutput("res"),
      dataTableOutput("similarDrugs")
    )
    )
))

