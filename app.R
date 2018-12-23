#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GoldenMutagenesis)
library(seqinr)

sequence_check<-function(input_sequence){
  input_sequence<-str_to_upper(input_sequence)
  if(nchar(input_sequence)%%3!=0) {
    stop(paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "))
  }
  codon_seq<-splitseq(s2c(input_sequence))
  met<-which(str_detect(codon_seq, "ATG"))
  if(length(met) == 0) {
    stop("No Methionine in the provided sequence. Stopping here. Please check the provided sequence.")
  }
  
  if(min(met) != 1){
    warning(paste("No Methionine at first codon found! Please check the provided sequence! Took codon #", min(met), "as start.", sep=" "))
    codon_seq<-codon_seq[min(met):length(codon_seq)]
  } #else(codon_seq<-codon_seq[-1])
  
  stop<-which(str_detect(codon_seq, "(TAA)|(TGA)|(TAG)"))
  if(length(stop) == 0) {
    stop("No stop codon in the provided sequence. Stopping here. Please check the provided sequence!")
  }
  
  if(max(stop) != length(codon_seq)) {
    warning(paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "))  
    codon_seq<-codon_seq[1:max(stop)]
  }# else {
  #codon_seq <- codon_seq[-length(codon_seq)]
  #}
  return(codon_seq)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  navbarPage("GoldenMutagenesis", id="Nav",
   tabPanel("Input and Settings",
   sidebarLayout(position = "right",
                 sidebarPanel(h1("Settings"), 
                              fluidRow(
                                column(10,
                                numericInput("temp", 
                                             h3("Target Temperature in Celsius"), 
                                             value = 60),
                                textInput("re", h3("Restriction Enzyme"), value = "GAAGAC")
                                )
                 )
                              ),
                 mainPanel(h1("Input"),
                           fluidRow(
                             column(3,
                                    textAreaInput("sequence", h3("Sequence"), width = '80%', resize = 'both', height='80%',
                                       value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA"),
                                    actionButton('next1', 'Domesticate')
                            )
                          )
                )
   )
   ), tabPanel("Mutations", 
               sidebarLayout(position = "right", 
                             sidebarPanel(h1("Selected Modifications"),
                                          fluidRow(
                                            column(3,h5("Position"), htmlOutput("position")),
                                            column(3,h5("Aminoacid") ,htmlOutput("aa"))
                                          )
                            ),
                            mainPanel(
                              h1("Selection"),
                              fluidRow(h3("Amino Acid Sequence"),
                                column(10, verbatimTextOutput("aa_sequence"))
                              ),
                              fluidRow(h3("Codon Sequence"), column(10, verbatimTextOutput("codon_sequence"))),
                              uiOutput("newmutation"),
                              fluidRow(column(3, actionButton('next2', 'Calculate Primers')
                                              )
                                       )
                            )
               )
      ),
   tabPanel("Primers", h3("Primer Report"), htmlOutput("primers"))
  )
)

#mutations<-c()
server <- function(input, output, session) {
  values <- reactiveValues()
  values$mutations<-c()
  observeEvent(input$next1, {
    updateTabsetPanel(session, "Nav",
                      selected = "Mutations")
  })
  observeEvent(input$next1, {
    values$mutations<-domesticate(input$sequence,restriction_enzyme = input$re);
    output$position<-renderUI({HTML(paste(sapply(values$mutations, function(x){return(x[1])}), collapse="<br>"))})
    output$aa<-renderUI({HTML(paste(sapply(values$mutations, function(x){return(x[2])}), collapse="<br>"))})
  })
  output$codon_sequence<-reactive(sequence_check(input$sequence))
  output$aa_sequence<-reactive(translate(s2c(input$sequence)))
  #output$position<-reactive(paste(sapply(mutations, function(x){return(x[1])}), collapse="<br>"))
  #output$aa<-reactive(paste(sapply(mutations, function(x){return(x[2])}), collapse="<br>"))
  output$newmutation<- renderUI({
    aminoacids<-a()
    names(aminoacids)<-aaa()
    fluidRow(
    column(2, selectInput("mpos", "Aminoacid Position", choices = 1:length(translate(s2c(input$sequence))))),
    column(2, selectInput("maa", "Aminoacid", choices = aminoacids), actionButton("mnew", label = "Add"))
    )
  })
  
  observeEvent(input$mnew, {
    values$mutations<-c(list(c(input$mpos, input$maa)), (values$mutations));
    output$position<-renderUI({HTML(paste(sapply(values$mutations, function(x){return(x[1])}), collapse="<br>"))})
    output$aa<-renderUI({HTML(paste(sapply(values$mutations, function(x){return(x[2])}), collapse="<br>"))})
  })
  
  observeEvent(input$next2, {
    updateTabsetPanel(session, "Nav",
                      selected = "Primers");
    print(values$mutations);
    primers<-mutate(input$sequence, replacements = values$mutations, restriction_enzyme = input$re, target_temp = input$temp)
    output$primers<-renderText(paste(capture.output(print_primer(primers)), collapse = "<br>"))
  })
}

# Run the application 


shinyApp(ui = ui, server = server)

