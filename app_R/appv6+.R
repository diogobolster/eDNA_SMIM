library(shiny)
library(shinyWidgets)
library(tidyverse)
library(DescTools)
# library(mathjaxr)
source("aux_fun.R")
library(optimx)
library(gt)
library(broom)
library(shinyhelper)
library(ggforce)
library(markdown)
library(ggnewscale)
library(expint)
library(shinyscreenshot)
# library(shinydashboard)

rm(list=ls())

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(div(
    column(6,strong("eDNA analysis tool")),
    column(6,tags$img(src="logo3.png",height="15%", width="15%", align = "right")),
  )),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(

      # helpText(h4("About:")),
      checkboxInput("about", "About", TRUE),
      # checkboxInput("diagram", "Ilustration", FALSE),
      helpText(h4("Data type:"))%>%
        helper(type = "markdown",
               content = "help_typedata"),
      h5("eDNA Concnetration"),
      radioButtons("data_type",label = NULL,
                   choices = list("total" = 1, "by size" = 2),
                   selected = 1,inline = T),
      
      helpText(h4("Input Data:")),
      checkboxInput("tab_data1", "Sampling Data", TRUE)%>% 
        helper(type = "markdown",
               content = "tab_data1"),
      # checkboxInput("tab_data2", "Sampling Data", FALSE)%>%
      #   helper(type = "markdown",
      #          content = "help"),
      
      uiOutput("checkbox_sourcedata"),
      

    helpText(h4("Analysis:")),
    uiOutput("LOD_analysis"),
    uiOutput("SMIM_analysis"),
    uiOutput("pred_analysis"),
    uiOutput("sens_analysis"),
    
    # checkboxInput("tab3", "Prediction", TRUE)%>% 
    #   helper(type = "markdown",
    #          content = "tab3help"),
    # checkboxInput("tab2", "Sensitivity", FALSE)%>% 
    #   helper(type = "markdown",
    #          content = "tab2help"),
    
    # helpText(h4("Total C Analysis:")),
    

    
    # checkboxInput("tab_totalC", "LOD", TRUE)%>% 
    #   helper(type = "markdown",
    #          content = "tab3help"),
    # checkboxInput("tab_SMIM", "SMIM", TRUE)%>% 
    #   helper(type = "markdown",
    #          content = "tab2help"),
    

    ,width = 2),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # textOutput('text')
      tabsetPanel(id = "tabs",
      
      ## tab_data2
      tabPanel("About",
         # br(), 
         # fluidRow(
         #   column(11,
         #  div(includeMarkdown("About.md"),style="text-align: justify;"))),
         # br(), 
         fluidRow(
         column(7,
                br(),
                div(includeMarkdown("About.md"),style="text-align: justify;"),
                # helpText(h4("Required information:")),
                # div(includeMarkdown("reqinf.md"),style="text-align: justify;"),
                # br(),
                # helpText(h5("How to Use:")),
                div(includeMarkdown("hw2use.md"),style="text-align: justify;"),
                # br(),
                # helpText(h4("Assumptions:")),
                includeMarkdown("assumptions.md"),
                # br(),
                # helpText(h4("Comments:")),
                div(includeMarkdown("comments.md"),style="text-align: justify;"),
                ),#end of column
         
         column(5,
                align="left",
                helpText(h4("Schematics:")),
                fluidRow(column(12,align="center",
                imageOutput("scheme", height = "190px"),
                h5("(a)"),
                imageOutput("scheme_SMIM", height = "190px"),
                h5("(b)"),
                imageOutput("scheme2", height = "190px"),
                h5("(c)"),
                imageOutput("scheme3", height = "240px"),
                h5("(d)"),
                ))
         ),#end of column
         ),#end of fluidrow


      ),            

# sampling data -----------------------------------------------------------

                  
       ## tab_data1   
      tabPanel("Sampling Data",
          
      fluidRow(
        # column(2,
               # radioButtons("data_type",label = "Select type of data:",
               #              choices = list("Total concentration" = 1, "Concentration per size" = 2),
               #              selected = 1)%>%
               #   helper(type = "markdown",
               #          content = "help#inputdata") ),
        column(2,
                 radioButtons("source",label = "Select type of input:",
                              choices = list("Upload Table" = 1, "Manual Input" = 2),
                              selected = 1)%>%
                   helper(type = "markdown",
                          content = "help#inputselection"),
        ),
        column(1,
               conditionalPanel(
                 condition = "input.source == 1",
               uiOutput("data_in"),
               )
        )
        
        ),
      
         conditionalPanel(
           condition = "input.source == 1",
           
           fluidRow(
             
           column(5,
             fluidRow(column(9,
                fluidRow(column(9,
              
              # Input select for file upload
              uiOutput("file_input"),
              # tags$style(HTML(".error-text { color: red; }")),
              textOutput("status")%>% 
                tagAppendAttributes( style = "color: red;")
             ),
             
             column(2,style='margin-left:-20px;', 
             downloadButton("downloadTemplate", "Template",
                            style = "font-size: 10px;"))
             ),
             uiOutput("treat_var_out"),
             
             # selectInput("treat_var", 
             #             label = h4("Treatment:"),
             #             choices = NULL)%>% 
             #   helper(type = "markdown",
             #          content = "treathelp"), 
             )),#end column # end of Fluid Row     
             uiOutput("regression_table"),
           
           # tableOutput("whatever"),
           # tableOutput(outputId = "regret")
           ),
           column(7,
           plotOutput(outputId = "regre")
           )), #end of column #endo of fluidrow
           # tableOutput(outputId = "raw")
         ),
         
        conditionalPanel(
           condition = "input.source == 2 &&  input.data_type == 2",
           fluidRow(
             column(3,textInput("treat1", label = "Treatment:",value='no_name')),
             column(2,
                    numericInput("numclass1", 
                                 label = "Number of classes:", value = 3)%>% 
                      helper(type = "markdown",
                             content = "help#classes"))
           ),
           fluidRow(
             column(5,
                    fluidRow(
                     column(7,
                            helpText(h4("Removal Rates:")%>% 
                                       helper(type = "markdown",
                                              content = "help#rates")),
                            fluidRow(
                              column(3,h5("Filter:")),
                              column(3,h5("Avg:")),
                              column(3,h5("LowerB:")),
                              column(3,h5("UpperB:")),
                            ),
                            uiOutput("dynamicUI")
                            ),
                     column(5,
                            helpText(h4("Intercept:"))%>%
                            helper(type = "markdown",
                                   content = "help#Intercept"),
                            fluidRow(
                              # column(3,h5("Filter:")),
                              column(4,h5("Avg:")),
                              column(4,h5("LowerB:")),
                              column(4,h5("UpperB:")),
                            ),
                            uiOutput("dynamicUI_Intercept")
                     )),
                    helpText(h4("Regression + 95% CI:")%>% 
                               helper(type = "markdown",
                                      content = "help#reg_table_man")),
                    
                    tableOutput("contents_man")),
              column(7,offset = 0,style='margin-left:0px;', 
                    plotOutput(outputId = "plot_man", height = "320px")
             ),# end of column
           ),# end of fluidrow
           
           
           
         ),

       ),  

# Source data -------------------------------------------------------------

       
      ## tab_data3
      tabPanel("Source Data",
           fluidRow(
             column(4,
                    # helpText(h4("Tipical Initial PSD:")),
                    fluidRow(column(6,
                    radioButtons("source_f0",label = "Select type of input:",
                                 choices = list("Database" = 1,
                                                "Upload source table" = 4,
                                                "Manual Input" = 2
                                                ),
                                 selected = 1)%>%
                      helper(type = "markdown",
                             content = "help#inputselection2"),
                    )),#end of fluidrow end of column
                    
              conditionalPanel(
                condition = "input.source_f0==1",
                # uiOutput("select_target"),
                selectInput("source_species_db",
                            label = "Species:",
                            choices = NULL),
                # uiOutput("slider_inputs_var")%>% 
                #   helper(type = "markdown",
                #          content = "dbhelp_source")
                # tableOutput('db')
                # textOutput("error_text")
              ),                    
                    
              conditionalPanel(
                condition = "input.source_f0==2",
                
                uiOutput("slider_inputs_var")%>% 
                  helper(type = "markdown",
                         content = "manhelp_source"),
                tableOutput('text_fcomp'),
                textOutput("error_text")
              ),

              conditionalPanel(
                condition = "input.source_f0==4",
                fluidRow(column(9,
                  fluidRow(column(9,                                
                    fileInput("file2", "Choose CSV File",
                              accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))%>% 
                      helper(type = "markdown",
                             content = "inputhelp_source")),
                    column(2,style='margin-left:-20px;', 
                           downloadButton("downloadTemplate2", "Template",
                                          style = "font-size: 10px;"))
                    ),
                
                  uiOutput("filter_source"),
                # radioButtons("filter_source",h4("Filter by:"),
                #              choices = list("Species" = 1, 
                #                             "Treatment" = 2
                #                             ),
                #              selected = 1)%>%
                #   helper(type = "markdown",
                #          content = "treathelp_source"),

                conditionalPanel(condition = "input.filter_source==1",
                  selectInput("source_species_var",
                              label = NULL,
                              choices = NULL)),
                
                
                
                conditionalPanel(condition = "input.filter_source==2",
                selectInput("source_treat_var",
                            label = NULL,
                            choices = NULL)))),
                
                
                
                
                
                
                # tableOutput("whatever_fspecies"),
              ),
              uiOutput("table_source"),
              # tableOutput("whatever_f1st")%>% 
              #   helper(type = "markdown",
              #          content = "tablehelp_source"),
              textOutput("error_text2"),
              
                    
             ),
             
             # column(1,tableOutput('text_fcomp')),
             column(8,
                    # conditionalPanel(condition = "input.add_sample > 0",
                    #                  radioButtons("ss2", h4("Sample:"),
                    #                               choices = list("Sample 1" = 1, "Sample 2" = 2),
                    #                               selected = 1, inline = T)),
                    helpText(h4("Particle size distribution:")),
                    # 
                    plotOutput("histPlot3"),heigh = 350)
           ),
               
      ),

# MC panel ----------------------------------------------------------------

   ##Tab 3
   # conditionalPanel(condition= "input.tab3==true",
   tabPanel("Prediction",
          
          
            fluidRow(
              column(3,
               div(
                 class = "box",
                 style = "border: 2px solid #000000; padding: 5px; margin-top: 5px;",      
                     (h4("Simulation set-up:"))%>% 
                   helper(type = "markdown", content = "MChelp"),
                     # conditionalPanel(condition = "input.add_sample > 0",
                     #                  radioButtons("ss3", h4("Sample:"),
                     #                               choices = list("Sample 1" = 1, "Sample 2" = 2),
                     #                               selected = 1, inline = T)),
                 actionButton(inputId = "submit",label = "Submit",
                              style="color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; width:230px;height:35px;font-size:140%;
                              font-weight:bold;text-align: center;"),
                 
                     numericInput('nsim', 'Number of extrapolation trials:', 1000, min = 1, max = 1e5),
                     numericInput('tmax', 'Maximum distance upstream (m):', 100, min = 30, max = 200),
                     numericInput('dt', 'Step(m):', 2, min = 30, max = 200),
                     h5("Tolerance (on source PSD):")%>% 
                   helper(type = "markdown", content = "tol_sourcehelp"),
                     numericInput('tol', label=NULL, 0.05, min = 0.001, max = 0.99),
                 h5("Cutoff filters:")%>% 
                   helper(type = "markdown", content = "help#cutoff"),
                 fluidRow(column(6,
                   div(numericInput('d_max', 'Distance(m):', NA, min = 0, max = 100)),
                   style = "font-size: 90%;"),
                   column(6,
                          # div(tableOutput("whateverMC"))
                 div(numericInput('c0_max', 'Released C0(cp/ml):', NA, min = 0, max = 100),
                    style = "font-size: 78%;")
                 )),
               # numericInput('c_bio', 'C/biomass factor:', NA, min = 0, max = 1e4)%>%
               #   helper(type = "markdown", content = "help#biom_fact")
               ),
               
                     
            ),
              column(9,
                     
                     
                     fluidRow(column(4,
                                    helpText(h4("Prediction distributions:"))%>% 
                                    helper(type = "markdown", content = "help#predic_dist"),
                                     h5(textOutput('nsolutions')),
                     )),
                     fluidRow(
                       column(3,
                              h5("Distance between source and 1st sample:")%>% 
                                helper(type = "markdown", content = "help"),
                              
                       ),
                       
                       column(3,
                              verbatimTextOutput('text3'),
                              tags$style("#text3 { font-size: 12px; font-weight: bold; }")
                              # numericInput('u_stream', 'U_average(m/hr):', 0, min = 0, max = 100),  
                              
                       ),
                       
                       column(3,
                              h5("Concentration released at source:")%>% 
                                helper(type = "markdown", content = "help"),
                              
                       ),
                       column(3,
                              # br(),
                              strong(verbatimTextOutput('text4')),
                              tags$style("#text4 { font-size: 12px; font-weight: bold; }")
                              
                       ),
                       
                     ),
                     
                     fluidRow(
                       column(6,
                              plotOutput(outputId = "tab3plothist",height = 350),
                              uiOutput("proxy1"),
                        # fluidRow(
                        #   column(2,
                        #          align="right",
                        #          h5("Proxy:")
                        #          ),
                        #   column(8,
                        #          align="left",
                        #          radioButtons("pdf", label = NULL,
                        #                       choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                        #                       selected = 1,inline = T)
                        #          )
                        # ),
                       ),
                       
                       column(6,
                              plotOutput(outputId = "tab3plothist2",height = 350),
                              uiOutput("proxy2"),
                              
                              # fluidRow(
                              #   column(2,
                              #          align="right",
                              #          h5("Proxy:")
                              #   ),
                              #   column(8,
                              #          align="left",
                              #          radioButtons("pdf2", label = NULL,
                              #                       choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                              #                       selected = 1,inline = T)
                              #   )
                              # ),
                              
                              
                              # radioButtons("pdf2", h5("pdf proxy:"),
                              #              choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                              #              selected = 1, inline = T)%>% 
                              #   helper(type = "markdown", content = "help")
                       ),
                       # column(2,radioButtons("pdf", h3("pdf proxy:"),
                       #                       choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                       #                       selected = 1))
                     ),        

            )),
            
            

            
            br(),
            fluidRow(
              column(4,
                     helpText(h4("Solution cases table:"))%>% 
                       helper(type = "markdown", content = "help#sol_table"),
                     div(tableOutput("whateverMC"),style = "font-size: 80%;"),
                     # div(tableOutput("whateverMC2"),style = "font-size: 65%;"),
                     ),
              column(8,
                     
                     helpText(h4("Back-extrapolations:"))%>% 
                       helper(type = "markdown", content = "help#sol_table"),
                     plotOutput("back_reg_mult",height = 350),
                     uiOutput("plot_MC",height = 350),
              )
              
            ),
            
            
                        
            fluidRow(
              column(6,uiOutput("plot_ui"))
            ),

            
   # ),
   ),


# single prediction -------------------------------------------------------


      ## tab 2
      # uiOutput("tab_panel2"),
      # conditionalPanel(condition= "input.tab2==true",
      tabPanel("Sensitivity",
               # conditionalPanel(condition= "input.tab2==true",
               fluidRow(
                 column(5,
                        div(
                          class = "box",
                          style = "border: 2px solid #000000; padding: 10px; margin-top: 10px;",       
                          
                          sliderInput(inputId = "time",
                                      label = h5("Distance|time before 1st sample (hr|m):"),
                                      min = 0,
                                      max = 100,
                                      value = 30)%>% 
                            helper(type = "markdown", content = "help"),
                          # tableOutput("whatever2"),
                          fluidRow(
                            column(6,
                                   withMathJax(
                                     helpText(h5("Removal rate range:")))%>% 
                                     helper(type = "markdown",
                                            content = "rem_rangehelp"),  #\\(\\lambda\\)
                                   uiOutput("slider_rem_inputs")),#end of column
                            
                            column(6,
                                   withMathJax(
                                     helpText(h5("Intercept range:")))%>% 
                                     helper(type = "markdown",
                                            content = "int_rangehelp"),
                                   uiOutput("slider_int_inputs")),#end of column
                            
                          )), # end of fluidrow
                        helpText(h5("Selection Table:")),
                        div(tableOutput('table_tab2'), style = "font-size: 90%; width: 100%"),
                        # helpText(h5("Start and End of back extrapolation:")),
                        # tableOutput("startendbe"),
                        
                 ), # end of column
                 column(7,
                        helpText(h4("Backward extrapolation:")),
                        plotOutput(outputId = "back_reg",height = 250),
                        # conditionalPanel(condition = "input.add_sample > 0",
                        # radioButtons("ss", h4("Sample:"),
                        #              choices = list("Sample 1" = 1, "Sample 2" = 2),
                        #              selected = 1, inline = T)),
                        helpText(h4("Particle size distribution:")),
                        plotOutput(outputId = "distPlot"),height = "250px")
               ),
               helpText(h4("Evolution in time|distance per size:")),
               fluidRow(
                 column(6,plotOutput(outputId = "normC")),
                 column(6,plotOutput(outputId = "frac"))
                 # column(6,plotOutput(outputId = "distPlot"))
               ),
               helpText(h4("Phase Diagram:")),
               fluidRow(
                 # column(6,plotOutput(outputId = "distPlot")),
                 column(6,plotOutput(outputId = "phase"))
               ),
               
               # ),
      ),## close tabPanel



# total C -----------------------------------------------------------------

tabPanel("LOD",
         
       # fluidRow(
       #   column(2,
       #          radioButtons("source_totalC",label = "Select type of input:",
       #                       choices = list("Upload Table" = 1, "Manual Input" = 2),
       #                       selected = 1)%>%
       #            helper(type = "markdown",
       #                   content = "help#inputselection_totalc"),
       #   )),
       
       
       # conditionalPanel(
       #   condition = "input.source_totalC == 1",
       #   
         fluidRow(

           column(4,
                  # fluidRow(column(9,
                    # fluidRow(column(9,
                    #                 fileInput("file_totalC", "Choose CSV File",
                    #                           accept = c(
                    #                             "text/csv",
                    #                             "text/comma-separated-values,text/plain",
                    #                             ".csv")
                    #                 )%>%
                    #                   helper(type = "markdown",
                    #                          content = "inputhelp_totalC")
                    # ),
       #   # 
                    # column(2,style='margin-left:-20px;',
                    #        downloadButton("downloadTemplate_totalc", "Template",
                    #                       style = "font-size: 10px;"))
                    # )),
       #   # 
                    # selectInput("treat_var_totalC",
                    #             label = h4("Treatment:"),
                    #             choices = NULL)%>%
                    #   helper(type = "markdown",
                    #          content = "treathelp_totalC"),
           # )),#end column # end of Fluid Row

                # helpText(h4("Regression + 95% CI:")%>%
                #            helper(type = "markdown",
                #                   content = "help#reg_table_totalC")),
                # tableOutput("contents_totalC")%>%
                #   helper(type = "markdown",
                #          content = "decayhelp"),
       
       # tableOutput("whatever"),
              
              numericInput('LOD', 'Limit of detection (LOD)(copies/ml)', 0.05, min = 1e-5, max = 1e5)%>%
           helper(type = "markdown",content = "lodhelp"),
              h5("Maximum detectable distance downstream (m):")%>%
              helper(type = "markdown", content = "lodhelp2"),
                tableOutput("raw_totalc"),
            # tableOutput("pred_totalc")
                # tableOutput(outputId = "regret")
         ),
           column(8,
                  # plotOutput(outputId = "regre_totalC"),
                  plotOutput(outputId = "plot_pred_totalC"),
                  
           )), #end of column #endo of fluidrow
         # tableOutput(outputId = "raw")
       # ),
       # # 
       # conditionalPanel(
       #   condition = "input.source_totalC == 2",
       #   fluidRow(
       #     column(3,textInput("treat1_totalC", label = "Treatment:",value='no_name')),
       #     # column(2,
       #     #        numericInput("numclass1",
       #     #                     label = "Number of classes:", value = 3)%>%
       #     #          helper(type = "markdown",
       #     #                 content = "help#classes"))
       #   ),
       #   fluidRow(
       #     column(5,
       #            fluidRow(
       #              column(7,
       #                     helpText(h4("Removal Rates:")%>%
       #                                helper(type = "markdown",
       #                                       content = "help#rates")),
       #                     fluidRow(
       #                       column(3,h5("Filter:")),
       #                       column(3,h5("Avg:")),
       #                       column(3,h5("LowerB:")),
       #                       column(3,h5("UpperB:")),
       #                     ),
       #                     # uiOutput("dynamicUI")
       #              ),
       #              column(5,
       #                     helpText(h4("Intercept:"))%>%
       #                       helper(type = "markdown",
       #                              content = "help#Intercept"),
       #                     fluidRow(
       #                       # column(3,h5("Filter:")),
       #                       column(4,h5("Avg:")),
       #                       column(4,h5("LowerB:")),
       #                       column(4,h5("UpperB:")),
       #                     ),
       #                     # uiOutput("dynamicUI_Intercept")
       #              )),
       #            helpText(h4("Regression + 95% CI:")%>%
       #                       helper(type = "markdown",
       #                              content = "help#reg_table_man")),
       # 
       #            tableOutput("contents_man_totalC")),
       #     column(7,offset = 0,style='margin-left:0px;',
       #            plotOutput(outputId = "plot_man_totalC", height = "320px")
       #     ),# end of column
       #   ),# end of fluidrow
       # 
       # 
       # 
       # ),
       # 
       
      ),## close tabPanel   


# SMIM --------------------------------------------------------------------

      tabPanel("SMIM",
                       
       fluidRow(column(3,
       conditionalPanel(condition = "input.data_type==2",                          
                       fluidRow(column(4,align="right",h5("Filter size:")),
                                column(5,align="left",
                                  selectInput("size2smim",
                                    label = NULL,
                                    choices = NULL))%>% 
                                  helper(type = "markdown", content = "help_SMIM_filter"),
       ),style='margin-bottom:-10px'),                       
         # numericInput('C0', label = h5("Denormalization factor (C0):"), 1, min = 1e-2, max = 1e6)%>%
         #   helper(type = "markdown", content = "helpx"),
         helpText(h4("Model parameters:"))%>% 
           helper(type = "markdown", content = "help_SMIM"),
         
         h5("Conservative parameters:")%>% 
           helper(type = "markdown", content = "help_SMIM_conservative"),
         
         fluidRow(
           column(3,
                  withMathJax(HTML('<h5><b>$$U(m/s):$$</b></h5>'))),
           # column(5,style="margin-top: 10px",
           #        htmlOutput('U'),
           #        # numericInput('U', label = NULL, 0.0982, min = 1, max = 1e6)
           # ),#end of column
           column(5,
                  numericInput('Us', label = NULL, 0.0982, min = 1, max = 1e6)),#end of column
         ),# end of fluidrow 
         
         fluidRow(
           column(3,
                  withMathJax(HTML('<h5><b>$$D(m2/s):$$</b></h5>'))),
           # column(5,style="margin-top: 10px",
           #        htmlOutput('D'),
           #        # numericInput('D', label = NULL, 0.0345, min = 1, max = 1e6)
           # ),#end of column
           column(5,
                  numericInput('Ds', label = NULL, 0.0345, min = 1, max = 1e6)),#end of column
         ),# end of fluidrow 
         
         fluidRow(
           column(3,
                  withMathJax(HTML('<h5><b>\\(\\lambda\\):</b></h5>'))
           ),
           # column(5,style="margin-top: 10px",
           #        htmlOutput('lambda'),
           #        # numericInput('lambda', label = NULL, 0.0850, min = 1, max = 1e6)
           # ),#end of column
           column(5,
                  numericInput('lambdas', label = NULL, 0.0850, min = 1, max = 1e6)),#end of column
         ),# end of fluidrow 
         
         fluidRow(
           column(3,
                  withMathJax(HTML('<h5><b>\\(\\beta\\):</b></h5>'))
           ),
           # column(5,style="margin-top: 10px",
           #        htmlOutput('beta'),
           #        # numericInput('beta', label = NULL, 0.8980, min = 1, max = 1e6)
           # ),#end of column
           column(5,
                  numericInput('betas', label = NULL, 0.8980, min = 1, max = 1e6)),#end of column
         ),# end of fluidrow 
         
         
         h5("Non-Conservative parameters:")%>% 
           helper(type = "markdown", content = "help_SMIM_nonconservative"),        
         
         fluidRow(
           column(3,
                  withMathJax(HTML('<h5><b>$$r_s(1/s):$$</b></h5>'))),
           column(5,style="margin-left: 0px",
                  htmlOutput('rs'),
                  numericInput('rss', label = NULL, 3e-5, min = 1, max = 1e6)
           ),#end of column
         ),#end of row
         
               fluidRow(
                 column(3,
                        helpText(h5(" "))),
                 column(4,
                        helpText(h5("Optimized:"))),#end of column
                 column(5,
                        helpText(h5("Manual:"))),#end of column
               ),# end of fluidrow  
               
               # fluidRow(
               #   column(3,
               #          HTML('<h5><b>rs(1/s):</b></h5>')),
               #   # column(4,style="margin-left: 10px",
               #   #        htmlOutput('rs'),
               #   #        numericInput('rs', label = NULL, 1e-3, min = 1, max = 1e6)
               #   #        ),#end of column
               #   column(5,
               #          numericInput('rss', label = NULL, 1e-3, min = 1, max = 1e6)),#end of column
               # ),# end of fluidrow 
               
               fluidRow(
                 column(3,
                        withMathJax(HTML('<h5><b>$$r_h(1/s):$$</b></h5>'))),
                 column(4,style="margin-top: 10px",
                        htmlOutput('rh'),
                        # numericInput('rh', label = NULL, 1e-3, min = 1, max = 1e6)
                        ),#end of column
                 column(5,
                        numericInput('rhs', label = NULL, 1e-3, min = 1, max = 1e6)),#end of column
               ),# end of fluidrow 
               
               
               # fluidRow(
               #   column(2,
               #          HTML('<h5><b>t1(s):10^ </b></h5>')),
               #   column(5,style="margin-top: 10px",
               #          htmlOutput('logT1'),
               #          # numericInput('logT1', label = NULL, 0, min = 1, max = 1e6)
               #          ),#end of column
               #   column(5,
               #          numericInput('logT1s', label = NULL, 0, min = 1, max = 1e6)),#end of column
               # ),# end of fluidrow 
               # 
               # fluidRow(
               #   column(2,
               #          HTML('<h5><b>t2(s):10^</b></h5>')),
               #   column(5,style="margin-top: 10px",
               #          htmlOutput('logT2'),
               #          # numericInput('logT2', label = NULL, 5, min = 1, max = 1e6)
               #          ),#end of column
               #   column(5,
               #          numericInput('logT2s', label = NULL, 5, min = 1, max = 1e6)),#end of column
               # ),# end of fluidrow 
            h5("Spatial Decline rate (1/m):")%>% 
           helper(type = "markdown", content = "help_declinerate"),
               fluidRow(
                 column(3,
                        withMathJax(HTML('<h5><b>$$k:$$</b></h5>'))),
                 column(4,style="margin-top: 10px",
                        htmlOutput('koptim')),#end of column
                 column(5,style="margin-top: 10px",
                        htmlOutput('ks')),#end of column
               ),# end of fluidrow 
         # tableOutput("whatever_optim"),
         
       ),#end of column
       column(9,
              fluidRow(
                column(3,
                       align = "left",
                       helpText(h4("Concentration profile:"))%>% 
                helper(type = "markdown", content = "help_SMIMprofile")
                )),
              fluidRow(
                column(2,
                       align = "right",
                       h5("Maximum distance downstream(m):")
                ),
                column(4,
                       sliderInput(inputId = "xlim",
                                   label=NULL,
                                   min = 0,
                                   max = 500,
                                   value = 100)
                ),
                column(3,
                         align = "right",
                         h5("Concentration Preditction at distance(m):")),
                  # HTML('<h5><b>xloc:</b></h5>')),
                  column(2,
                         numericInput('xloc', label = NULL, 50, min = 1, max = 1e6)
                         ),#end of column
                  # column(4,
                  #        numericInput('xlocs',label = NULL, 100, min = 1, max = 1e6)),#end of column
                
                
                ),
                plotOutput("SMIM_Plot"),
              
              
              # fluidRow(
              #   column(4,
              #          verbatimTextOutput('text_smim_1'),
              #          tags$style("#text_smim_1 { font-size: 14px; font-weight: bold; }")
              #   ),
              #   
              #   
              #   column(4,
              #          strong(verbatimTextOutput('text_smim_2')),
              #          tags$style("#text_smim_2 { font-size: 14px; font-weight: bold; }")
              #   ),
              # ),
              # h4("SMIM model schematic:", a(href="https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JG004684", "(Roche 2019)")),
              # HTML("This is some help text with a <a href='https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13751' target='_blank'>Roche 2019</a>."),
              # imageOutput("img")
              )
       )# end of fluidrow         
                        
      ),## close tabPanel  
      ),## close tabsetPanel

      
    width = 10) ## close main panel
    ), ## close sidebar layout
setBackgroundColor("#d3d3d3")
)## close fluid page




server <- function(input, output,session) {
  
  observe_helpers(help_dir = "help_mds",withMathJax = TRUE)
  

  # observeEvent(input$about, {
  #   screenshot(scale=3)
  # })

# about_tab ---------------------------------------------------------------

  output$scheme <- renderImage({
    
    list(src = "LOD_schematic.png",
         width = "90%",
         height = 200)
    
  }, deleteFile = F)
  
  output$scheme_SMIM <- renderImage({
    
    list(src = "SMIM_schematic2.png",
         contentType = 'image/png',
         width = "90%",
         height = 200)
    
  }, deleteFile = F)
  
  output$scheme2 <- renderImage({
    
    list(src = "schematic2s2.png",
         contentType = 'image/png',
         width = "90%",
         height = 200)
    
  }, deleteFile = F)
  
  output$scheme3 <- renderImage({
    
    list(src = "plot_eDNA_components3.png",
         contentType = 'image/png',
         width = "80%",
         height = 250)
    
  }, deleteFile = F)
  
# data_tabs ---------------------------------------------------------------
  # control of appearance of the tab
  observeEvent(input$about, {
    if(input$about==FALSE){
      hideTab(inputId = "tabs", target = "About")
    }
    else{
      showTab(inputId = "tabs", target = "About")
    }
  })
  
  
  observeEvent(input$tab_data1, {
    if(input$tab_data1==FALSE){
      hideTab(inputId = "tabs", target = "Sampling Data")
    }
    else{
      showTab(inputId = "tabs", target = "Sampling Data")
    }
  })

  observeEvent(input$tab_data2, {
    if(input$tab_data2==FALSE){
      hideTab(inputId = "tabs", target = "Sampling Data")
    }
    else{
      showTab(inputId = "tabs", target = "Sampling Data")
    }
  })
  
  
  observeEvent(input$tab_data3, {
    if(input$tab_data3==FALSE){
      hideTab(inputId = "tabs", target = "Source Data")
    }
    else{
      showTab(inputId = "tabs", target = "Source Data")
    }
  })
  
  
  observeEvent(input$tab_totalC, {
    if(input$tab_totalC==FALSE){
      hideTab(inputId = "tabs", target = "LOD")
    }
    else{
      showTab(inputId = "tabs", target = "LOD")  
    }
  })
  
  observeEvent(input$tab_SMIM, {
    if(input$tab_SMIM==FALSE){
      hideTab(inputId = "tabs", target = "SMIM")
    }
    else{
      showTab(inputId = "tabs", target = "SMIM")  
    }
  })
  
  
  output$downloadTemplate <- downloadHandler(
    filename = function() {
      "template_sampling.csv"
    },
    content = function(file) {
      # Copy the template file to the download location
      if(input$data_type==2){
        if(input$data_in==1){
        
      file.copy("template_sampling_perSize.csv", file)
        }
        else{
      file.copy("template_sampling_perSize_hours.csv", file)
        }
        }
      if(input$data_type==1){
        if(input$data_in==1){
        file.copy("template_sampling_Total.csv", file)
        }
        else{
        file.copy("template_sampling_Total_hours.csv", file)  
        }
        }
      
    }
  )
  
  output$downloadTemplate2 <- downloadHandler(
    filename = function() {
      "template_source.csv"
    },
    content = function(file) {
      # Copy the template file to the download location
      file.copy("template_source.csv", file)
    }
  )

# check csv format table --------------------------------------------------

  observeEvent(input$file1, {
    if (!is.null(input$file1)) {
      df <- read.csv(input$file1$datapath)
      if(input$data_type==2){ 
        if(input$data_in==1){
      expected_columns <- 
        c("Treatment", "Replicate", "Filter","Target","Concentration","Distance") 
      expected_num_columns <- length(expected_columns)
        }
        else{
          expected_columns <- 
            c("Treatment", "Replicate", "Filter","Target","Concentration","Hours") 
          expected_num_columns <- length(expected_columns)
        }
      }
      else{
        if(input$data_in==1){
        expected_columns <- 
          c("Treatment", "Replicate","Target","Concentration","Distance")
        expected_num_columns <- length(expected_columns)
        }
        else{
        expected_columns <- 
          c("Treatment", "Replicate","Target","Concentration","Hours")
        expected_num_columns <- length(expected_columns)
        }
      }
      
      if (!all(colnames(df) == expected_columns)) {
        output$status <- renderText({
          paste("Incorrect column names. Expected:", paste(expected_columns, collapse = ", "))
        })
        return()
      }
      
      # Check if number of columns matches
      if (ncol(df) != expected_num_columns) {
        output$status <- renderText({
          paste("Incorrect number of columns. Expected:", expected_num_columns, "columns.")
        })
        return()
      }
      
      uploaded_data(read.csv(input$file1$datapath))
      val(10)
      output$status <- renderText({
        " "
      })
      
      
      
    }})
  
# reading/regression tab --------------------------------------------------

  
  # Dynamic UI for fileInput
  
  output$data_in <- renderUI({
    radioButtons("data_in",label = "Data in:",
                 choices = list("Space" = 1, "Time" = 2),
                 selected = 1)%>%
      helper(type = "markdown",
             content = "help#datain")
  })
  
  output$file_input <- renderUI({
    fileInput("file1", "Upload CSV File",accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv"))
  })
  
  
  # Reactive value to store uploaded data
  uploaded_data <- reactiveVal(NULL)
  val <- reactiveVal(NULL)
  val2 <- reactiveVal(NULL)
  
  totalC <- reactiveVal(1)
  perSize <- reactiveVal(NULL)
  # Event handler for file upload and input change
  intime<- reactiveVal(NULL)
  
  
  observeEvent(input$file1, {
    if (!is.null(input$file1)) {
      # Read the uploaded CSV file
      req(val())
      
      output$treat_var_out <- renderUI({
        selectInput("treat_var", 
                    label = h4("Treatment:"),
                    choices = NULL)%>% 
          helper(type = "markdown",
                 content = "treathelp")
      })
      
      output$regression_table <- renderUI({
        helpText(h4("Regression + 95% CI:")%>% 
                   helper(type = "markdown",
                          content = "help#reg_table"))
        
        tableOutput("contents")%>% 
          helper(type = "markdown",
                 content = "decayhelp")
      })
      # output$checkbox_sourcedata <- renderUI({
      # checkboxInput("tab_data3", "Source Data", FALSE)%>%
      #   helper(type = "markdown",
      #          content = "tab_data2")
      # })
    }
  })
  
  observeEvent(input$data_in,{
    uploaded_data(NULL)
    val2(NULL)
    
    output$file_input <- renderUI({
      fileInput("file1", "Upload CSV File")
    })
    
    output$treat_var_out <- renderUI({
    })
    
    output$regression_table <- renderUI({
    })
    
    if(input$data_in==2){
      intime(10)
      hideTab(inputId = "tabs", target = "SMIM")
      output$SMIM_analysis <- renderUI({  
        # checkboxInput("tab_SMIM", "SMIM", FALSE)%>% 
        #   helper(type = "markdown",
        #          content = "tab2help")
      })
     }
    else{
      
      showTab(inputId = "tabs", target = "SMIM")
      output$SMIM_analysis <- renderUI({  
        checkboxInput("tab_SMIM", "SMIM-rates", TRUE)%>% 
          helper(type = "markdown",
                 content = "tab2help")
        })
      
      }
    
    
  })

  
  observeEvent(input$data_type, {
    # Reset uploaded data to NULL
    uploaded_data(NULL)
    val(NULL)
    val2(NULL)
    # data_regre(NULL)
    # Remove the existing fileInput and render a new one
    output$file_input <- renderUI({
      fileInput("file1", "Upload CSV File")
    })
    
    output$status <- renderText({
      " "
    })
    
    output$treat_var_out <- renderUI({
    })
    
    output$regression_table <- renderUI({
    })
    
    if(input$data_type==2){
      
      totalC(NULL)
      perSize(1)
      
    showTab(inputId = "tabs", target = "Source Data")
    output$checkbox_sourcedata <- renderUI({
      checkboxInput("tab_data3", "Source Data", TRUE)%>% 
        helper(type = "markdown",
               content = "tab_data2")
    })
    showTab(inputId = "tabs", target = "Prediction")
    output$pred_analysis <- renderUI({
    checkboxInput("tab3", "Prediction", TRUE)%>% 
      helper(type = "markdown",
             content = "tab3help")
    })
    
    hideTab(inputId = "tabs", target = "Sensitivity")
    output$sens_analysis <- renderUI({
    checkboxInput("tab2", "Sensitivity", FALSE)%>% 
      helper(type = "markdown",
             content = "tab2help")
    })
    
    # hideTab(inputId = "tabs", target = "LOD")
    # output$LOD_analysis <- renderUI({
    # # checkboxInput("tab_totalC", "LOD", FALSE)%>% 
    # #   helper(type = "markdown",
    # #          content = "tab3help")
    # })  
    # hideTab(inputId = "tabs", target = "SMIM")
    #   output$SMIM_analysis <- renderUI({  
    # # checkboxInput("tab_SMIM", "SMIM", FALSE)%>% 
    # #   helper(type = "markdown",
    # #          content = "tab2help")
    #   })
    }
    
    if(input$data_type==1){
      totalC(1)
      perSize(NULL)
      
      hideTab(inputId = "tabs", target = "Source Data")
      output$checkbox_sourcedata <- renderUI({
        # checkboxInput("tab_data3", "Source Data", FALSE)%>% 
        #   helper(type = "markdown",
        #          content = "tab_data2")
      })
      hideTab(inputId = "tabs", target = "Prediction")
      output$pred_analysis <- renderUI({
        
        # checkboxInput("tab3", "Prediction", FALSE)%>% 
        #   helper(type = "markdown",
        #          content = "tab3help")
      })
      
      hideTab(inputId = "tabs", target = "Sensitivity")
      output$sens_analysis <- renderUI({
        # checkboxInput("tab2", "Sensitivity", FALSE)%>% 
        #   helper(type = "markdown",
        #          content = "tab2help")
      })
      
      showTab(inputId = "tabs", target = "LOD")
      output$LOD_analysis <- renderUI({
        checkboxInput("tab_totalC", "LOD-distance", TRUE)%>% 
          helper(type = "markdown",
                 content = "LOD_help")
      })  
      
      if(is.null(intime())){
      showTab(inputId = "tabs", target = "SMIM")
      output$SMIM_analysis <- renderUI({  
        checkboxInput("tab_SMIM", "SMIM-rates", TRUE)%>% 
          helper(type = "markdown",
                 content = "SMIM_help")
      })
      }
    }
    
    
  })
  
  observeEvent(input$submit,{
  output$proxy1 <- renderUI({
  fluidRow(  
  column(2,
         align="right",
         h5("Proxy:")
  ),
  column(8,
         align="left",
         radioButtons("pdf", label = NULL,
                      choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                      selected = 3,inline = T)
  )
  )
  })
  })
  
  observeEvent(input$submit,{
  output$proxy2 <- renderUI({  
  fluidRow(
    column(2,
           align="right",
           h5("Proxy:")
    ),
    column(8,
           align="left",
           radioButtons("pdf2", label = NULL,
                        choices = list("No" = 1, "Normal" = 2, "Lognormal" = 3),
                        selected = 1,inline = T)
    )
  )
  })
  })    
  
  # reading csv file
  dat <- reactive({
    req(input$file1)
    req(uploaded_data())
    
    # uploaded_data(read.csv(input$file1$datapath))
    data <- uploaded_data() #read.csv(input$file1$datapath)
    names <- unique(data$Treatment)
    
    if(input$data_in==2){
      
      data = data%>%rename(Distance=Hours)
    }
    
    if(input$data_type==2){
    sizes <- sort(as.numeric(unique(data$Filter)))
    num_parameters <- length(unique(data$Filter))

    list(data = data, names = names, sizes=sizes, num_parameters = num_parameters)
    }
    else{
      list(data = data, names = names)
    }
})
  
  # Choosing treatment
  observeEvent(input$source,{
    if(input$source==2){
      updateSelectInput(session, "treat_var", choices = input$treat1)
    }
    else{
      observeEvent(dat(), {
        updateSelectInput(session, "treat_var", choices = c(dat()$names,"ALL"))
        
      })
    }
    
  })
  # observeEvent(val(), {
  # # req(val())
  # 
  # })
  # 

  # Manual input dynamic interface
  output$dynamicUI <- renderUI({
      num <- input$numclass1
      if(num<2){
        stop("minimum of 2 classes")
      }
      else{
      # default values for manual input  
      fval=c(0.2,0.8,3,10,30,100)
      aval=c(0.1,0.1,0.1,0.25,0.3,0.05)
      lval=aval-0.07
      hval=aval+0.07
      
      inputs <- lapply(1:num, function(i) {
        fluidRow(
        column(3,offset = 0,style='padding:1px;',
               numericInput(paste0("mfilter", i), 
                     NULL, value = fval[i])),         #paste0("Filter ", i, ":")
        column(3,offset = 0,style='padding:1px;'
               ,numericInput(paste0("mavg", i), 
                              NULL, value = aval[i])),   #paste0("Avg", i, ":")
        column(3,offset = 0,style='padding:1px;',
               numericInput(paste0("mlow", i), 
                     NULL, value = lval[i])),      # paste0("Low", i, ":")
        column(3,offset = 0,style='padding:1px;',
               numericInput(paste0("mhigh", i), 
                     NULL, value = hval[i]))    #paste0("High", i, ":")
        )
      })
      }
      do.call(tagList, inputs)
  })
  # getting Intercept manual inputs
  output$dynamicUI_Intercept <- renderUI({
    num <- input$numclass1
    if(num<2){
      stop("minimum of 2 classes")
    }
    else{
      # default values for manual input  
      # fval=c(0.2,0.8,3,10,30,100)
      aval=-c(0,0,0,0.25,0.3,0.05)
      lval=aval-0.3
      hval=aval+0.3
      
      inputs <- lapply(1:num, function(i) {
        fluidRow(
          # column(3,numericInput(paste0("mfilter", i), 
          #                       paste0("Filter ", i, ":"), value = fval[i])),
          column(4,offset = 0,style='padding-right:1px;',
                 numericInput(paste0("iavg", i), 
                              NULL, value = aval[i])),
          column(4,offset = 0,style='padding:1px;',
                 numericInput(paste0("ilow", i), 
                              NULL, value = lval[i])),
          column(4,offset = 0,style='padding:1px;',
                 numericInput(paste0("ihigh", i), 
                                NULL, value = hval[i]))
        )
      })
    }
    do.call(tagList, inputs)
  })

  
  # Manual input saved on a dataframe format
  mdat <-reactive({
    req(input$mfilter1)
    Nsizes=input$numclass1
    Treatment= input$treat1
    Filter=as.numeric(lapply(1:(Nsizes), function(i) {
      input[[paste0("mfilter", i)]]}))
    estimate = as.numeric(lapply(1:(Nsizes), function(i) {
      -input[[paste0("mavg", i)]]}))
    term="Distance"
    conf.low=as.numeric(lapply(1:(Nsizes), function(i) {
      -input[[paste0("mhigh", i)]]}))
    conf.high=as.numeric(lapply(1:(Nsizes), function(i) {
      -input[[paste0("mlow", i)]]}))

    a = tibble(Treatment,Filter,estimate,term,conf.low,conf.high)

    Filter=as.numeric(lapply(1:(Nsizes), function(i) {
      input[[paste0("mfilter", i)]]}))
    estimate = as.numeric(lapply(1:(Nsizes), function(i) {
      input[[paste0("iavg", i)]]}))
    term="(Intercept)"
    conf.high=as.numeric(lapply(1:(Nsizes), function(i) {
      input[[paste0("ihigh", i)]]}))
    conf.low=as.numeric(lapply(1:(Nsizes), function(i) {
      input[[paste0("ilow", i)]]}))

    b =rbind(a,tibble(Treatment,Filter,term,estimate,conf.low,conf.high))%>%
      arrange(Filter,term)

  })
  

  ## number of size classes
  Nsizes<- reactive({
    Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
  })
  
  
  ## preparing data for regressions
  raw.data.n<-reactive({
    req(input$file1)
    req(val())
    
    raw.data = dat()$data
    raw.data = filter(raw.data, Distance>=0)
    
    if(input$data_type==2){
    raw.data$Filter  = as.factor(raw.data$Filter)
    raw.data$Treatment  = as.factor(raw.data$Treatment)
    raw.data$Replicate  = as.factor(raw.data$Replicate)
    raw.data$Distance  = as.numeric(raw.data$Distance)
    
    ## normalizing mesocosm data with max value
    raw.data.n=raw.data %>% group_by(Treatment,Replicate) %>%
      mutate(max_c = max(Concentration,na.rm = TRUE))%>%
      mutate(norm_c = Concentration/max(Concentration,na.rm = TRUE)) %>%
      mutate(ln_c = ifelse((Concentration==0),NA,log(Concentration)))%>%
      mutate(ln_norm_c = ifelse((norm_c==0),NA,log(norm_c)))
    ## calculating the fractions based on normlaized values
    raw.data.n = raw.data.n %>% group_by(Treatment,Replicate,Distance)%>%
      mutate(frac = norm_c/sum(norm_c)) %>%ungroup()
    
    raw.data.n
    }
    else{
    
    
    # raw.data$Filter  = as.factor(raw.data$Filter)
    raw.data$Treatment  = as.factor(raw.data$Treatment)
    raw.data$Replicate  = as.factor(raw.data$Replicate)
    raw.data$Distance  = as.numeric(raw.data$Distance)
    
    ## normalizing mesocosm data with max value
    raw.data.n=raw.data %>% group_by(Treatment,Replicate) %>%
      mutate(max_c = max(Concentration,na.rm = TRUE))%>%
      mutate(norm_c = Concentration/max(Concentration,na.rm = TRUE)) %>%
      mutate(ln_c = ifelse((Concentration==0),NA,log(Concentration)))%>%
      mutate(ln_norm_c = ifelse((norm_c==0),NA,log(norm_c)))
    ## calculating the fractions based on normlaized values
    # raw.data.n = raw.data.n %>% group_by(Treatment,Replicate,Distance)%>%
    #   mutate(frac = norm_c/sum(norm_c)) %>%ungroup()
    
    raw.data.n 
    }
    })
  # 
  denorm_factor<- reactive({
    req(input$source)
    if(input$source==1){
    req(input$treat_var)
    a = raw.data.n()%>%subset(Treatment==input$treat_var)%>%
      summarise(average_maxC=mean(max_c,na.rm = TRUE))
    a$average_maxC
    }
    else if(input$source==2){
      1
    }
  })
  output$raw <- renderTable({
    as.numeric(denorm_factor())*c(1,1,0.1)
    # raw.data.n()%>%subset(Treatment==input$treat_var)%>%
    #   mutate(average_maxC=mean(max_c,na.rm = TRUE))##%>%
    # subset(Filter==10)
  })
  
  ## regressions
  table_var<-reactive({
    req(input$treat_var)
    r = raw.data.n()
    if(input$treat_var=='ALL'){
      r$Treatment="ALL"
      # nested <- raw.data.n() %>% nest (data = -c(Filter)) %>%
      #   arrange(Filter)
    }
    else{ }
    
    if(input$data_type==2){
      nested <- r %>% nest (data = -c(Treatment,Filter)) %>%
        arrange(Treatment,Filter)
    }
    else{
      nested <- r %>% nest (data = -c(Treatment)) %>%
        arrange(Treatment)  
    }
    # }
    
    regressions <-nested %>% 
      mutate(
        fit = map(data, possibly(~ lm(ln_norm_c~Distance,data=.x),NA)),
        tidied = map(fit,conf.int = TRUE, tidy),
      )
    
    
    if(input$data_type==2){
    table = regressions %>% unnest(tidied) %>% #filter(term=="Distance") %>%
      select(Treatment,Filter,term, estimate,std.error,p.value,conf.low,conf.high) %>%
      {if(input$treat_var!='ALL') subset(.,Treatment==input$treat_var) else .}
    }
    else{
      table = regressions %>% unnest(tidied) %>% #filter(term=="Distance") %>%
        select(Treatment,term, estimate,std.error,p.value,conf.low,conf.high) %>%
        {if(input$treat_var!='ALL') subset(.,Treatment==input$treat_var) else .}  
    }
    
  })
  
  ## output table with avearhe estimates and 95% confidence interval
  output$contents <- renderTable({
    req(val())
    t = table_var()%>%ungroup()%>%
      mutate(Avg = estimate,.after = term)%>%
      select(-c(Treatment,p.value,estimate,std.error))
    t$term[t$term=='Distance']='removal rate'
    t
    # table_var()
      
      # mutate(Low = -conf.high,
      #        High = -conf.low)%>%select(-c(estimate,std.error,conf.low,conf.high))
  })

  output$contents_man <- renderTable({
    # req(val())
    t = mdat()
    t$term[t$term=='Distance']='removal rate'
    t%>%select(-Treatment)
    # mdat()%>%select(Filter,estimate,term)%>%
    #   pivot_wider(values_from = estimate,names_from = term)%>%
    #   rename(Intercept="(Intercept)",slope=Distance)
    # mdat()
  })
  
  # filtering data for only one treatment
  data_regre<- eventReactive(input$treat_var,{
    raw.data.n() %>% #subset(Treatment==input$treat_var)
      {if(input$treat_var!='ALL') subset(.,Treatment==input$treat_var) else .}
  })
  
  observeEvent(data_regre(),{
    val2(10)
  })
  
  # plot of regressions for one selected treatment
  output$regre <- renderPlot({
    # req(input$file1)
    req(val2())
    
      # subset(Filter==10)%>%
      gg <- ggplot(data_regre())+
      geom_point(aes(x=Distance, y=ln_norm_c),size=2.5,shape=21)+
      geom_smooth(aes(x=Distance, y=ln_norm_c),method='lm')+
      # if(){facet_grid(cols = vars(Filter))+}
      # else{}
      theme_bw()+
      theme(text = element_text(size = 15))+
      ylab('ln(normalized(copies/ml))')+
      labs(x='Distance (m)')
      
      if (input$data_type==2) {
        gg <- gg + facet_grid(cols = vars(Filter))
      }
      
      if(input$data_in==2){
        gg<-gg + labs(x='Time (hr)')
      }
      gg
  })  
  
  # plot of regressions for one selected treatment
  output$plot_man <- renderPlot({
    req(input$numclass1)
    t=mdat()
    t$estimate[t$term=='Distance']=-mdat()$estimate[mdat()$term=='Distance']
    t$conf.low[t$term=='Distance']=-mdat()$conf.high[mdat()$term=='Distance']
    t$conf.high[t$term=='Distance']=-mdat()$conf.low[mdat()$term=='Distance']
    
    
    avg = t%>%select(Filter,estimate,term)%>%
      pivot_wider(values_from = estimate,names_from = term)%>%
      rename(Intercept="(Intercept)",slope=Distance)
    
    low = t%>%select(Filter,conf.low,term)%>%
      pivot_wider(values_from = conf.low,names_from = term)%>%
      rename(Intercept="(Intercept)",slope=Distance)
    
    
    high = t%>%select(Filter,conf.high,term)%>%
      pivot_wider(values_from = conf.high,names_from = term)%>%
      rename(Intercept="(Intercept)",slope=Distance)
    
      ggplot()+
      geom_abline(data=avg,aes(slope=-slope, intercept=Intercept,color='Avg'),size=2)+
      geom_abline(data=low,aes(slope=-slope, intercept=Intercept,color='LB'),size=2)+
      geom_abline(data=high,aes(slope=-slope, intercept=Intercept,color='UB'),size=2)+
      # geom_smooth(aes(x=Distance, y=ln_norm_c),method='lm')+
      facet_grid(cols = vars(Filter))+
      theme_bw()+
      theme(text = element_text(size = 15),legend.position = c(0.87, 0.80))+
      ylab('ln(copies/ml)')+
      scale_y_continuous(limits=c(-1,1))+
      scale_x_continuous(expand=c(0,0),limits=c(0,10))+
      # expand_limits(x = 0, y = -1)+
        # xlim(0,10)+
      labs(x='Distance (hr | m)')+
      scale_color_manual(name=NULL,
        values = c("LB" = "red","UB" = "blue","Avg"="black")
      )  
  })  
  
  

# source database ---------------------------------------------------------------- 
db <- reactive({
  Filter <- c(0.2,1,10)
  STHD <- c(0.31,0.49,0.21)
  CARP <- c(0.22,0.16,0.62)
  
  db = tibble(Filter,STHD,CARP)%>%
    pivot_longer(-c(Filter),names_to='target',values_to = "mean")%>%
    mutate(std=0.0001)
})
  
  output$db<- renderTable({
    pts()
  })

  observeEvent(db(), {
    updateSelectInput(session, "source_species_db", 
                      choices = unique(db()$target))})
  
  

# source data -------------------------------------------------------------
  # read file input for source PSD
  source_dat <- reactive({
    req(input$file2)
    data <- read.csv(input$file2$datapath)
    names <- unique(data$Treatment)
    species <- unique(data$Target)
    sizes <- sort(as.numeric(unique(data$Filter)))
    num_parameters <- length(unique(data$Filter))
  #   
     list(data = data, names = names, species=species, sizes=sizes, num_parameters = num_parameters)
  })
  
  observeEvent(input$file2, {
    if (!is.null(input$file2)) {
  
  output$filter_source <- renderUI({
    radioButtons("filter_source",h4("Filter by:"),
                 choices = list("Species" = 1, 
                                "Treatment" = 2
                 ),
                 selected = 1)%>%
      helper(type = "markdown",
             content = "treathelp_source")
  })
  }
  })
  
  observeEvent(pts(), {
    req(perSize())
    output$table_source <-renderUI({
      tableOutput("whatever_f1st")%>% 
        helper(type = "markdown",
               content = "tablehelp_source")
    })
  })

  observeEvent(source_dat(), {
    updateSelectInput(session, "source_treat_var", choices = source_dat()$names)})
  
  observeEvent(source_dat(), {
    updateSelectInput(session, "source_species_var", choices = source_dat()$species)})


  
  data_source <-reactive({
    req(input$filter_source)
    if(input$filter_source==1){
      source_dat()$data%>%subset(Target==input$source_species_var)
      }
    else{
      source_dat()$data%>%filter(Treatment==input$source_treat_var)
    }
    
  })
  
  fr0_fromsourcetable<- reactive({
    source_dat()$data%>%filter(Treatment==input$source_treat_var)%>%
      mutate(fr0=Concentration/sum(Concentration))%>%arrange(Filter)
  })
  
  fr0_species<- reactive({
    req(input$source_f0)
    # source_dat()$data%>%subset(Target==input$source_species_var)%>%
    
    data_source()%>%
      group_by(Treatment,Replicate)%>%mutate(frac=Concentration/sum(Concentration))%>%
      ungroup()%>%group_by(Filter)%>%
      summarise(mean=mean(frac,na.rm=TRUE),std=sd(frac,na.rm=TRUE))  
    
  })

  #manual input
  output$slider_inputs_var <- renderUI({
    sliderInputs <- lapply(1:(Nsizes()-1), function(i) {
      # paste0(),
      numericInput(paste0("fr0", i), paste0("Filter: ", filters()[i]), min = 0, max = 1, value = 1/Nsizes())
    })
    tags$div(
      style = "display: flex; flex-direction: column; justify-content: center;",
      do.call(tagList, sliderInputs)
    )
  })
  #manual input complementary
  fman <-reactive({
    fr=as.numeric(lapply(1:(Nsizes()-1), function(i) {
      input[[paste0("fr0", i)]]+1e-5}))
    fman=tibble(Filter=filters(),mean=c(fr,1-sum(fr)))%>%mutate(std=0.0001)
    # fman$std=0.001                 
  })
  
  #manual input complementary
  output$text_fcomp <- renderTable({
    df=tibble(fman()$mean[Nsizes()])
    names(df) = paste0("Filter: ", filters()[Nsizes()])
    df
  })
  
  #manual input error
  output$error_text <- renderText({
    # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
    req(input$fr01)
    if (fman()$mean[Nsizes()]<0) {
      stop("fractions should add to 1")
    } else {
    }
  })
  
  
  # final values for initial PSD 
  f_initial <- reactive({
    req(input$source_f0)
    req(perSize())
    
    if(input$source_f0==2){
      req(input$fr01)
      ft=fman()}
    # else if(input$source_f0==1){
    #   ft=fr0_fromtable()$fr0}
    else if(input$source_f0==4){
      req(input$source_species_var)
      ft=fr0_species()}
    else{
      req(input$source_species_db)
      if(input$data_type==2){
      ft = db()%>%subset(target==input$source_species_db)%>%
      select(-target)}}
    # else{ft = as.numeric(cinput()$f)}
    # return(ft)
  }) 
  
  f_1st_sample <-reactive({
    req(perSize())
    if(input$source==1){
    req(val2())  
    data_regre()%>%subset(Distance==0 & frac<1)%>%
      group_by(Filter)%>%summarise(mean=mean(frac,na.rm=TRUE),std=sd(frac,na.rm=TRUE))}
    else{
      mdat()%>%subset(term=="(Intercept)")%>%
        pivot_longer(c(estimate,conf.low,conf.high),
                     names_to='names',values_to = "ln_norm_c")%>%
        mutate(Distance=0,norm_c=exp(ln_norm_c))%>%
        group_by(names)%>%mutate(frac=norm_c/sum(norm_c))%>%
        ungroup()%>%group_by(Filter)%>%
        summarise(mean=mean(frac,na.rm=TRUE),std=sd(frac,na.rm=TRUE))
    }
  })
  
  pts<-reactive({
    req(f_1st_sample())
    req(f_initial())
    req(perSize())
    
    
    df1=f_initial()%>%mutate(row='source')
    # df1$row = 'source'
    
    df2=f_1st_sample()%>%mutate(row='1st_sample')
    # df1$row = 'sample'
    pts=rbind(df1,df2)
    # if(input$source==1){
    # psdt = f_1st_sample()$mean
    # else{psdt = f_1st_sample()$mean}
    # psdt= cselec2()/sum(cselec2())
    
    # pts = tibble()
    # pts <- pts%>%rbind(f_initial(),psdt)
    # 
    # for(i in 1:Nsizes()){
    #   colnames(pts)[i] = paste0("f", i)
    # }
    # pts <- pts%>%mutate(row=c('i','f'))
  })
  

  # output$whatever_fspecies <- renderTable({
  #   # f_initial()C
  #   # f_1st_sample()%>%rename(frac=mean)
  #   # fr0_species()
  #   pts()
  #   # df1=f_initial()%>%mutate(row='source')
  #   # df2=f_1st_sample()%>%mutate(row='sample')
  #   # pts=rbind(df1,df2)
  # })
    
  output$whatever_f1st <- renderTable({
    # source_dat()$data
    # f_1st_sample()
    # f_1st_sample()%>%rename(frac=mean)
    # fr0_species()
    
    
    pts()%>%rename(fraction=mean,error=std,location=row)
    # pts()%>%rename_with(~ filters(),-c(row))%>%
    #   pivot_longer(cols = -c(row), values_to = "f", names_to = "size")
  })
  
  output$error_text2 <- renderText({
    # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
    req
    
    if (sum(f_1st_sample()$mean)<0.99  | sum(f_1st_sample()$mean)>1.01) 
      {
      stop("fractions should add to 1 (Check Sampling Data Input)")
      } 
    else { }
  })
  

 
  
  output$histPlot3 <- renderPlot({
    
    pts()%>%#rename_with(~ filters(),-c(row))%>%
      #pivot_longer(cols = -c(row), values_to = "f", names_to = "size")%>%
      # ggplot(aes(x = factor(as.numeric(size),levels=filters()), y = f, fill = factor(row)))+
      ggplot(aes(x = factor(Filter), y = mean, fill = factor(row)))+
      geom_col(position=position_dodge2(padding = 0,reverse = TRUE))+
      geom_errorbar(aes(ymin = mean-std, ymax = mean+std),
                    position=position_dodge2(width = .3,padding = 0.2,reverse = TRUE))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      # ylim(0,1)+
      scale_fill_manual(name = "", 
                        breaks= c("source", "1st_sample"),
                        labels = c("Source", "1st_Sample"),
                        values=c("darkgreen",
                                 "orange"))+
      xlab('Filter')+ylab('fraction')
    
    
    
  })
  
  
  # cselec2 <- eventReactive(input$ss2,{
  #   if(input$ss2==1){cinput()$c}
  #   else{cinput2()$c2}
  # })
  # 
  # cselec3 <- eventReactive(input$ss3,{
  #   if(input$ss3==1){cinput()$c}
  #   else{cinput2()$c2}
  # })
  
  
  
  
# Single prediction -------------------------------------------------------------------

  # control of appearance of the tab
  observeEvent(input$tab2, {
    if(input$tab2==FALSE){
    hideTab(inputId = "tabs", target = "Sensitivity")
    }
    else{
    showTab(inputId = "tabs", target = "Sensitivity")  
    }
  })
  
  # choosing between manual vs csv input for sizes store filters names
  filters <- reactive({
    if(input$source==1) {
      filters=dat()$sizes #csv input
      }
    
    else{
      filters=unique(mdat()$Filter) # manual input
      }  
    
  })
  
  #removal rates multipliers
  output$slider_rem_inputs <- renderUI({
    
    sliderInputs <- lapply(1:Nsizes(), function(i) {
      sliderInput(paste0("mlbd_", i), paste0("filter:", filters()[i]), min = 0, max = 1, value = 0.5)
    })
    tags$div(
      style = "display: flex; flex-direction: column; justify-content: center;",
      do.call(tagList, sliderInputs)
    )
  })
  
  #Interception multipliers
  output$slider_int_inputs <- renderUI({
    
    sliderInputs <- lapply(1:Nsizes(), function(i) {
      sliderInput(paste0("mint_", i), paste0("filter:", filters()[i]), min = 0, max = 1, value = 0.5)
    })
    tags$div(
      style = "display: flex; flex-direction: column; justify-content: center;",
      do.call(tagList, sliderInputs)
    )
  })
  
  # may not be necessary
  # from manual sampling tab
  # cselec <- eventReactive(input$ss,{
  # if(input$ss==1){cinput()$c}
  #   else{cinput2()$c2}
  # })
  
  # inputs from sliders
  # inputs<-reactive({
  #   tibble(ct=cselec(),
  #   mlbd=lapply(1:Nsizes(), function(i) {
  #       1-input[[paste0("mlbd_", i)]]
  #     }),
  #   lbd=lapply(1:Nsizes(), function(i) {
  #     input[[paste0("lbd_", i)]]
  #   })
  #   )
  #   
  # })
  
  multi <-reactive({
    mlbd=lapply(1:Nsizes(), function(i) {
      1-input[[paste0("mlbd_", i)]]
    })
    mint=lapply(1:Nsizes(), function(i) {
      input[[paste0("mint_", i)]]
    })
    tibble(mint,mlbd)%>%pivot_longer(everything())
  })

  #interpolation of multipliers
  table_pred <-reactive({
    req(input$mlbd_1,input$mint_1)
    if(input$source==1){table=table_var()}
    else{table=mdat()
         }
    
    table%>%
    select(Treatment,Filter,term,estimate,conf.low,conf.high)%>%
    mutate(multi = multi()$value)%>% 
    mutate(.,lambda = apply(.,1,function(x) 
      as.numeric(unlist(x[5]))+(as.numeric(unlist(x[6]))-as.numeric(unlist(x[5])))
      *as.numeric(unlist(x[7]))))
    # 
    # mutate(
    #  psdt= ct/sum(ct),
    #  c0=apply(.,1,function(x) as.numeric(x[6])*exp(-as.numeric(x[9])*input$time)),
    #  psd0=c0/sum(c0))
  })
  
  #single prediction table
  pred <-reactive({
    
    lambda =  table_pred()$lambda
    # filt =  levels(table_pred()$Filter)
    # length(filt)
    tmax = input$time
    dt = 1
    
    t = single_pred(tmax,dt,lambda,filters())
    
  })
  output$table_tab2 <- renderTable({
    t = table_pred()%>%ungroup%>% 
      select(Filter,term,estimate,conf.low,conf.high,lambda)%>%
      rename(selection = lambda)
    t$term[t$term=='Distance']='removal rate'
    t
    # %>%select(-c(mlbd,lbd))%>%
    #   mutate(Avg = -estimate,.after = Filter)%>%
    #   mutate(Low = -conf.high,
    #          High = -conf.low,.after = Avg)%>%
    #   mutate(Rem_rate = -lambda,.after=High)%>%
    #   mutate(c0=formatC(c0, digits = 3))%>%
    #   select(-c(estimate,conf.low,conf.high,lambda))%>%
    #   rename(C_sample=ct,
    #          PSD_Sample=psdt,
    #          C_Predict=c0,
    #          PSD_Predict=psd0)
    # 
  })
  
  output$startendbe <- renderTable({
    
    # pred() %>% select(time,Filter,frac)%>%
    #   mutate(Filter=paste0("f_",as.factor(Filter)))%>%
    #   pivot_wider(names_from = Filter,values_from = frac)
    pred() %>% subset(time==0 | time == -input$time)%>%arrange(time)%>%
      rename('time|distance' = time)
    # mdat()%>%subset(term=="(Intercept)")%>%
    #   pivot_longer(c(estimate,conf.low,conf.high),
    #                names_to='names',values_to = "ln_norm_c")%>%
    #   mutate(Distance=0)
  })
  
  
  
  
  # pred_var <-reactive({
  #   Distance<-seq(0, input$time, length = 100)
  #   pred<-tibble(Distance=Distance)
  #   
  #   for(i in 1:Nsizes()){
  #     col_name = paste0("C", i)
  #     pred[[col_name]]<-df()$c0[i]*exp(df()$lambda[i]*Distance)
  #   }
  #   pred
  # })
    

  
  
  output$distPlot <- renderPlot({

    data = pred() %>% subset(time==0 | time == -input$time)
      # pivot_longer(cols = c(psdt,psd0), values_to = "psd", names_to = "time") %>%
      gg = ggplot(data,aes(x = factor(Filter), y = frac, fill = factor(time)))+
      geom_col(position = "dodge")+
      theme_bw()+
      theme(text = element_text(size = 15),legend.position = c(0.87, 0.80))+
      # ylim(0,1)+
      scale_fill_manual(name = "distance (m)", 
                        labels = c(paste0("@d=",-input$time),paste0("@d=0")),
                        values=c("darkgreen",
                                 "orange"))+
      xlab('Filter')+ylab('fraction')
    
    if(input$data_in==2){
      gg<-gg + 
        scale_fill_manual(name = "time (hr)", 
                          labels = c(paste0("@t=",-input$time),paste0("@t=0")),
                          values=c("darkgreen",
                                   "orange"))
        
    }
    gg
    
    
    
  },height= 300)
  
  

  
# 
#   frac_var <-reactive({
#     # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
#     df <- pred_var()%>%mutate(total= rowSums(select(.,-Distance)))
#     
#     for(i in 1:Nsizes()){
#       col_name <- paste0("f", i)
#       # df <- df %>% mutate(!!col_name := .[,i+1]/total)
#       df[[col_name]] = df[,i+1]/df$total
#       df[[col_name]] <- as.numeric(unlist(df[[col_name]]))
#     }
#     df
#   })
  
  output$back_reg <- renderPlot({
    
    if(input$source==1){data=data_regre()}
    else{data=mdat()%>%subset(term=="(Intercept)")%>%
      pivot_longer(c(estimate,conf.low,conf.high),
                   names_to='names',values_to = "ln_norm_c")%>%
      mutate(Distance=0)
      }
    
    
      # subset(Filter==10)%>%
      gg <- ggplot(data)+
      geom_point(aes(x=Distance, y=ln_norm_c),size=2.5,shape=21)+
      geom_smooth(aes(x=Distance, y=ln_norm_c),method='lm')+
      geom_line(data=pred(),aes(x=time, y=ln_norm_c))+
      facet_grid(cols = vars(Filter))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      ylab('ln(norm(copies/ml))')+
      labs(x='Distance (m)')
    
    
    if(input$data_in==2){
      gg<-gg + labs(x='Time (hr)')
    }
    gg
    
  })
  
  output$normC <- renderPlot({
    
    if(input$source==1){data=data_regre()}
    else{data=mdat()%>%subset(term=="(Intercept)")%>%
      pivot_longer(c(estimate,conf.low,conf.high),
                   names_to='names',values_to = "ln_norm_c")%>%
      mutate(Distance=0)%>%mutate(norm_c=exp(ln_norm_c))
    }
    # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
    my_colors <- rainbow(Nsizes())
    # pred()%>%
      # pivot_longer(cols = -c(Distance), values_to = "C", names_to = "Size")%>%
      gg = ggplot(pred(),aes(x = time, y = norm_c, color = factor(Filter)))+geom_line(size=2)+
      theme_bw()+
      # scale_y_continuous(trans='log10')+
      theme(text = element_text(size = 15),legend.position = c(0.87, 0.75))+
      scale_color_manual(name = "Filter",
                         labels =  unique(table_pred()$Filter),
                         values= my_colors)+
      labs(y="Norm_Concentration(copies/ml)")+
      geom_point(data=data,aes(x=Distance, y=norm_c, color = factor(Filter)),size=2.5,shape=21)+
      labs(x='Distance (m)')
    
    if(input$data_in==2){
      gg<-gg + labs(x='Time (hr)')
    }
    gg
    
  })
  

  
  
  output$frac <- renderPlot({
    # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
    my_colors <- rainbow(Nsizes())
    
      # pivot_longer(cols = -c(Distance), values_to = "C", names_to = "Size")%>%
      gg = ggplot(pred(),aes(x = time, y = frac, color = factor(Filter)))+geom_line(size=2)+
      theme_bw()+
      # scale_y_continuous(trans='log10')+
      theme(text = element_text(size = 15),legend.position = c(0.87, 0.75))+
      scale_color_manual(name = "Filter",
                         labels =  unique(table_pred()$Filter),
                         values= my_colors)+
      labs(y="Fraction")+
      {if(input$source==1)
        geom_point(data=data_regre(),aes(x=Distance, y=frac, color = Filter),size=2.5,shape=21)}+
      labs(x='Distance (m)')
    
    if(input$data_in==2){
      gg<-gg + labs(x='Time (hr)')
    }
    gg
  })

  # output$tab1plot3 <- renderPlot({
  #   # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
  #   my_colors <- rainbow(Nsizes())
  #   frac_var() %>%select(Distance,matches("^f"))%>%
  #     pivot_longer(cols = -c(Distance), values_to = "f", names_to = "Size")%>%
  #     ggplot(aes(x = Distance, y = f, color = Size))+geom_line(size=2)+
  #     theme_bw()+
  #     # scale_y_continuous(trans='log10')+
  #     theme(text = element_text(size = 15),legend.position = c(0.87, 0.75))+
  #     scale_color_manual(name = "Size",
  #                        labels =  df()$Filter,
  #                        values= my_colors)+
  #     labs(y="Fraction")
  # 
  # })
  # 
  #
  output$phase <- renderPlot({
    
    
    filt = paste0("f_",levels(as.factor(pred()$Filter)))
    
    
    df = pred() %>% select(time,Filter,frac)%>%
      mutate(Filter=paste0("f_",as.factor(Filter)))%>%
      pivot_wider(names_from = Filter,values_from = frac)%>%
      pivot_longer(cols=c(filt[2:Nsizes()]),values_to = 'f',names_to = 'Filter')
    
    my_colors <- rainbow(Nsizes())
    # df = pred() %>% select(time,Filter,frac)%>%
    #   mutate(Filter=paste0("f_",as.factor(Filter)))%>%
    #   pivot_wider(names_from = Filter,values_from = frac)
    
    
      gg = ggplot(df)+
      geom_point(aes_string(x = filt[1] , y = 'f', color = 'time', size = 'Filter'))+
      theme_bw()+
      # scale_y_continuous(trans='log10')+
      theme(text = element_text(size = 15))+
      geom_abline(intercept=1,slope=-1)+
      labs(y= "fn")+
        labs(color = "Distance \n (m)")
      
      if(input$data_in==2){
        gg<-gg + labs(color = "Time \n (hr)")
      }
      gg
      
      })
  
  
  
  # output$tab1plot4 <- renderPlot({
  #   # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
  #   my_colors <- rainbow(Nsizes())
  #   frac_var() %>%select(Distance,matches("^f"))%>%
  #     pivot_longer(cols = -c(Distance,f1), values_to = "f", names_to = "Size")%>%
  #     ggplot()+
  #     geom_point(aes(x = f1, y = f, color = Distance, size=Size))+
  #     theme_bw()+
  #     # scale_y_continuous(trans='log10')+
  #     theme(text = element_text(size = 15),legend.position = c(0.87, 0.7))+
  #     labs(x= "f1", y= "fn")
  #     # scale_fill_gradientn(colours = terrain.colors(10))
  #   # scale_color_manual(name = "Size",                       # labels   scale_color_manual(name = "Size",
  #   #                    labels =  df()$Filter[-(1)],
  #   #                    values= my_colors[-(1)])
  # })
  

# tabMC --------------------------------------------------------------------
  

  observeEvent(input$tab3, {
    if(input$tab3==FALSE){
      hideTab(inputId = "tabs", target = "Prediction")
    }
    else{
      showTab(inputId = "tabs", target = "Prediction")  
    }
  })
  
  # MC simulations
  dfMC <- reactive({
    req(perSize())
    if(input$source==1){table=table_var()}
    else{table=mdat()}
    
    table = table%>%
      select(Filter,term,estimate,conf.low,conf.high)
    
    withProgress(message = 'Calculation in progress', {
    MC_new(table,input$nsim,input$dt,input$tmax,Nsizes())
    })
    
  }) %>% bindEvent(input$submit)
  

  # back extrapolation plot
  output$back_reg_mult <- renderPlot({
    req(perSize())
    dfsol = dfMC()%>%
      filter(Nsim %in% df_close()$Nsim)%>%
      filter(Nsim<=1000)
    
    if(input$source==1){data=data_regre()%>%subset(Treatment==input$treat_var)}
    else{data=mdat()%>%subset(term=="(Intercept)")%>%
      pivot_longer(c(estimate,conf.low,conf.high),
                   names_to='names',values_to = "ln_norm_c")%>%
      mutate(Distance=0)
    }
    
    df = data%>% subset(Treatment==input$treat_var)
      # subset(Filter==10)%>%
      gg = ggplot(df)+
      geom_point(aes(x=Distance, y=ln_norm_c),size=2.5,shape=21)+
      geom_smooth(aes(x=Distance, y=ln_norm_c),method='lm')+
      geom_line(data=dfMC(),aes(x=time, y=ln_norm_c,color=as.factor(Nsim)))+
      scale_color_manual(values = rep("gray",input$nsim))+
      new_scale_colour()+
      geom_line(data=dfsol,aes(x=time, y=ln_norm_c,color=as.factor(Nsim)))+
      scale_color_manual(values = rep("blue",input$nsim))+
    # labs(color = "Possible Cases(time)")+
      facet_grid(cols = vars(Filter))+
      theme_bw()+
      theme(text = element_text(size = 15),legend.position = "none")+
      ylab('ln(norm(copies/ml))')+
      labs(x='distance (m)')
      
      if(input$data_in==2){
        gg<-gg + labs(x='Time (hr)')
      }
      gg
      

    
  })
  
  
  # Filtering the possible cases/solutions
  df_close <- reactive({
    req(perSize())
    dfMC()%>% select(time,Filter,frac,Nsim)%>%
      mutate(Filter=paste0("f_",as.factor(Filter)))%>%
      pivot_wider(names_from = Filter,values_from = frac)%>%
      mutate(SSe = apply(.,1,function(x) 
        sqrt(sum((x[3:(Nsizes()+1)]-f_initial()$mean[1:(Nsizes()-1)])^2))))%>%
      group_by(Nsim)%>%mutate(minSSe=min(SSe))%>%ungroup()%>% 
      filter(minSSe==SSe)%>%filter(SSe<=input$tol)%>% 
      filter(time!=-input$tmax,time!=0)%>%
      mutate(time_sim=paste(time,"_",Nsim))%>%
      {if(is.numeric(input$d_max) & input$d_max>0) filter(.,time>=-input$d_max) else .}
  }) 
  
  df_close_C0<- reactive({
    dfMC()%>%mutate(time_sim=paste(time,"_",Nsim))%>%
      filter(time_sim %in% df_close()$time_sim)%>%
      select(time,Filter,norm_c,Nsim)%>%
      mutate(Filter=paste0("f_",as.factor(Filter)))%>%
      pivot_wider(names_from = Filter,values_from = norm_c)%>%
      mutate(C0 = apply(.,1,function(x) sum(x[3:(Nsizes()+2)])))%>%
      mutate(factor = denorm_factor())%>%
      mutate(C0denorm = C0*factor)%>%
      {if(is.numeric(input$c0_max) & input$c0_max>0) filter(.,C0denorm<=input$c0_max) else .}
    
    # %>%
    #   {if(is.numeric(input$c_bio) & input$c_bio>0) mutate(.,bio=C0/input$c_bio) else .}
    # browser()
  })
  
  # Filtering the possible cases/solutions entire timeseries
  df_close_full <- reactive({
    req(perSize())
    dfMC()%>%
      filter(Nsim %in% df_close_C0()$Nsim)%>%
      select(time,Filter,frac,Nsim)%>%
      mutate(Filter=paste0("f_",as.factor(Filter)))%>%
      pivot_wider(names_from = Filter,values_from = frac)
      
  })
  
# Outputs   
  # Table of possible solutions
  output$whateverMC <- renderTable({
    df = df_close()%>%
    filter(Nsim %in% df_close_C0()$Nsim)%>%
    select(-c(time_sim,minSSe))%>%
    mutate(C0=formatC(as.integer(df_close_C0()$C0denorm),
      format="e",digits = 2),.after = time)%>%         
      slice_head(n=15)%>%rename('distance'="time")
      
    if(input$data_in==2){
      df%>%rename(time=distance)
    }
    else{
      df
    }
    # browser()
    # dfMC()
  })
  
  
  
  
  # output$whateverMC2 <- renderTable({
  #   df_C0()%>%slice_head(n=15)%>%select(c(time,Nsim,C0,factor,C0denorm))
  #   # %>%select(-c(time_sim,minSSe))%>%
  #   #   slice_head(n=15)%>%rename('distance'="time")
  #   
  #   # dfMC()
  # })
  # 
  output$nsolutions <- renderText({
    paste('Possible solutions: ', nrow(df_close_C0()))
  })
    # pts3<-reactive({
    #   # if(input$source_f0==2){ft=ft()}
    #   # else{ft=fr0_fromtable()$fr0}     
    #   # Nsizes=ifelse(input$source==1, dat()$num_parameters, input$numclass1)
    #   psdt= cselec3()/sum(cselec3())
    #   pts = tibble()
    #   pts <- pts%>%rbind(f_initial(),psdt)
    #   
    #   for(i in 1:Nsizes()){
    #     colnames(pts)[i] = paste0("f", i)
    #   }
    #   pts <- pts%>%mutate(row=c('i','f'))
    # })
    # 
    # 
    # output$tablepoints <- renderTable({
    #     pts3()%>%mutate(time=c('Initial','Sample'))%>%select(-row)
    # })
    
  output$plot_MC1 <- renderPlot({
      filt = paste0("f_",levels(as.factor(dfMC()$Filter)))
      
      pts = pts()%>%mutate(Filter=factor(paste0("f_",Filter)))%>%
        select(-std)%>%pivot_wider(names_from = Filter,values_from = mean)
      
      pts2 = pts%>%filter(row=='source')
      
      df = dfMC()%>%filter(Nsim<=1000)%>%
        select(time,Filter,frac,Nsim)%>%
        mutate(Filter=paste0("f_",as.factor(Filter)))%>%
        pivot_wider(names_from = Filter,values_from = frac)
      
        p = ggplot(df)+
        geom_point(aes_string(x=filt[1],y=filt[2],color='time'))+
        scale_color_gradient(high = "lightgray",low = "azure4") +
        labs(color = "All cases\n (Distance(m))")+
        
        geom_abline(intercept=1,slope=-1)+
        ylim(0,1)+xlim(0,1)+theme_bw()+
        theme(text = element_text(size = 15))+
        geom_point(data=pts,aes_string(x=filt[1],y=filt[2],fill = 'row'),pch = 21,size=3)+
        ylab(paste0('fraction of  ', filt[2], ' mm'))+
        xlab(paste0('fraction of  ', filt[1], ' mm'))+
        labs(title = "Phase diagram:")+
        scale_fill_manual(name='',
                          labels = c('Source PSD', '1st Sample PSD'),
                          values = c('source'='darkgreen',  '1st_sample'='orange'))+
        geom_circle(data=pts2,aes_string(x0=filt[1],y0=filt[2],r='input$tol'),
                    color = "red", fill = "lightgray", alpha = 0.2)
        
        if(input$data_in==2){
          p = p + labs(color = "All cases\n (Time(hr))")
        }
      
      if(nrow(df_close_full())>0){
        p = p + new_scale_colour()+
        geom_point(data=df_close_full(),aes_string(x=filt[1],y=filt[2],color='time'))+
        labs(color = "Possible Cases \n (Distance(m))")+
        geom_point(data=pts,aes_string(x=filt[1],y=filt[2],fill = 'row'),pch = 21,size=3)+
        geom_circle(data=pts2,aes_string(x0=filt[1],y0=filt[2],r='input$tol'),
                      color = "red", fill = "lightgray", alpha = 0.2)
        
        
        if(input$data_in==2){
          p = p + labs(color = "Possible Cases \n (Time(hr))")
        }

      }
        
        

      p
        
      
    })
    
  # plot organizer for more than 3 sizes    
  output$plot_MC <- renderUI({
      if(Nsizes()<=3){
      fluidRow(
        column(12, plotOutput("plot_MC1"))
      )
      }
      else{
        fluidRow(
          column(6, plotOutput("plot_MC1")),
          column(6, plotOutput("plot_MC2"))
        )
      }
    })%>% bindEvent(input$submit)
    
        
    # UI for plot based on condition
    # output$plot_ui <- renderUI({
    #   if(Nsizes()>4){
    #     plotOutput("tab3plot3")
    #   }
    # })%>% bindEvent(input$submit)
    # 
    # output$tab3plot3 <- renderPlot({
    # 
    #     if(input$u_stream>0){
    #       dfMC_d()%>%filter(Nsim<=1000)%>%
    #         ggplot()+geom_point(aes(x=f01,y=f04,color=-distance))+
    #         labs(color = "Distance")+
    #         geom_abline(intercept=1,slope=-1)+
    #         ylim(0,1)+xlim(0,1)+theme_bw()+
    #         theme(text = element_text(size = 15),legend.position = c(0.86, 0.75))+
    #         geom_point(data=pts3(),aes(x=f1,y=f4,fill = row),pch = 21,size=3)+
    #         ylab(paste0('f04:', "(Size:", filters()[4],")"))+
    #         xlab(paste0('f01:', "(Size:", filters()[1],")"))+
    #         labs(title = "Simulations:")+
    #         scale_fill_manual(name='',
    #                           labels = c('Source PSD', 'Sampled PSD'),
    #                           values = c('i'='darkgreen',  'f'='orange')
    #         )+
    #         geom_circle(data=pts3(),aes(x0=f1[1],y0=f4[1],r=input$tol),
    #                     color = "red", fill = "lightgray", alpha = 0.2) 
    #       
    #     }
    #     else{
    #       
    #       dfMC_d()%>%filter(Nsim<=1000)%>%
    #         ggplot()+geom_point(aes(x=f01,y=f04,color=-time))+
    #         labs(color = "Time")+
    #         geom_abline(intercept=1,slope=-1)+
    #         ylim(0,1)+xlim(0,1)+theme_bw()+
    #         theme(text = element_text(size = 15),legend.position = c(0.86, 0.75))+
    #         geom_point(data=pts3(),aes(x=f1,y=f4,fill = row),pch = 21,size=3)+
    #         ylab(paste0('f04:', "(Size:", filters()[4],")"))+
    #         xlab(paste0('f01:', "(Size:", filters()[1],")"))+
    #         labs(title = "Simulations:")+
    #         scale_fill_manual(name='',
    #                           labels = c('Source PSD', 'Sampled PSD'),
    #                           values = c('i'='darkgreen',  'f'='orange')
    #         )+
    #         geom_circle(data=pts3(),aes(x0=f1[1],y0=f4[1],r=input$tol),
    #                     color = "red", fill = "lightgray", alpha = 0.2)  
    #         
    #       
    #     }  
    # 
    # })%>% bindEvent(input$submit)
    
    output$text3 <- renderText({

      x = -df_close_C0()$time
      
      if(input$data_in==1){
      paste(#"Based on",length(x),"cases:\n",
            #"Range:",suppressWarnings(min(x)),"-",
            #suppressWarnings(max(x)),"hs|m\n",
            "Mean:",round(mean(x),2),"m\n",
            "SD:",round(sd(x),2),"m\n",
            "Mode:",ifelse(length(x)<=1,'NaN',round(mode(x),1)),"m"
            # "Mode:",round(Mode(x)[1],2)
      )
      }
      
      else{
        paste(#"Based on",length(x),"cases:\n",
          #"Range:",suppressWarnings(min(x)),"-",
          #suppressWarnings(max(x)),"hs|m\n",
          "Mean:",round(mean(x),2),"hr\n",
          "SD:",round(sd(x),2),"hr\n",
          "Mode:",ifelse(length(x)<=1,'NaN',round(mode(x),1)),"hr"
          # "Mode:",round(Mode(x)[1],2)
        ) 
      }
      
    }) 
    

    

    output$tab3plothist <- renderPlot({
      req(nrow(df_close_C0())>0)
      x = -df_close_C0()$time
      df=tibble(x)
       p =  ggplot() + 
          geom_histogram(data=df,aes(x=x,y=..density..),color="black",
                         breaks=seq(0,input$tmax,input$tmax/15), fill="darkgreen")+
          theme_bw()+theme(text = element_text(size = 15))+
          xlim(-input$tmax/20,input$tmax)+
         labs(x='Distance (m)')
       
       if(input$data_in==2){
         p = p + labs(x='Time (hr)')
       }
      
      if (input$pdf == 1){}
       
      else{
      if (input$pdf == 2){
        p = p + stat_function(fun=dnorm,
                              args=list(mean=mean(x),
                                        sd=sd(x)),
                              col='red',
                              size=2)
        
      }
      else{
        x = x + 1e0
        p = p + stat_function(fun=dlnorm,
                             args=list(meanlog=mean(log(x)),
                                       sdlog=sd(log(x))),
                             col='red',
                             size=2)
         
      }
      } 
       p
       
    })
    


    output$text4 <- renderText({
      # if(is.numeric(input$c_bio) & input$c_bio>0) {x = df_C0()$bio}
      # else {x = df_C0()$C0}
      x = df_close_C0()$C0denorm
      # if(input$source==1) {x = x*as.numeric(denorm_factor())}
      # else { }
      x=log10(x)
      paste(#"Based on",length(x),"cases:\n",
            #"Range:",suppressWarnings(formatC(min(x), format = "e", digits = 2)),"-",
            #suppressWarnings(formatC(max(x),format = "e", digits = 2)),"copies/ml | Kg\n",
            "Mean: 10^",formatC(mean(x), digits = 2),"cp/ml \n",
            "SD: 10^",formatC(sd(x), digits = 2),"cp/ml\n",
            "Mode: 10^",ifelse(length(x)<=1,'NaN',formatC(mode(x), digits = 2)),"cp/ml"
      )
      
    })
    
    output$tab3plothist2 <- renderPlot({
      req(nrow(df_close_C0())>0)
      # if(is.numeric(input$c_bio) & input$c_bio>0) {x = df_C0()$bio}
      # else {x = df_C0()$C0}
      x = df_close_C0()$C0denorm
      # if(input$source==1) {x = x*as.numeric(denorm_factor())}
      # else { }
        
        x=log10(x)
        df=tibble(x)
        x_seq <- seq(min(x)-1, max(x)+1, length.out = 100)
        y_seq <- dnorm(x_seq, mean = mean(x), sd = sd(x))
        y_seq2 <- dlnorm(x_seq, meanlog = mean(log(x)), sdlog = sd(log(x)))
        df2 = tibble(x_seq,y_seq,y_seq2)
        
        p =  ggplot(data=df) + 
          geom_histogram(aes(x=x,y=..density..),color="black",
                         bins=15,
          # breaks = seq(-bin_size, max(df_close()$Nfish), by = bin_size),
          # breaks=seq(mean(df_close()$C0)-3*sd(df_close()$C0),mean(df_close()$C0)+3*sd(df_close()$C0),sd(df_close()$C0)/5),
          fill="darkorange")+ #
          theme_bw()+theme(text = element_text(size = 15))+
          xlim(min(x)-1,max(x)+1)+
          labs(x='log10(Concentration (copies/ml))')

        if (input$pdf2 == 1){ }

        else{
          if (input$pdf2 == 2){
            p = p +
              geom_line(data = df2, aes(x = x_seq, y = y_seq),
                        color = "red", size=2)
          }
          else{
            p = p +
              geom_line(data = df2, aes(x = x_seq, y = y_seq2),
                        color = "red", size=2)
          }
        }
        p
    })
    

# tabrESULTSTABLE --------------------------------------------------------------------


    
    observeEvent(input$tab4, {
      if(input$tab4==FALSE){
        hideTab(inputId = "tabs", target = "Results Table")
      }
      else{
        showTab(inputId = "tabs", target = "Results Table")  
      }
    })
    
    output$table3 <- renderDataTable({
      df_close()%>%select(-c(minSSe))%>%mutate(across(everything(), round, digits=2))
    },options = list(paging = FALSE,searching = FALSE,dom='t'))
    
    
    output$downloadData <- downloadHandler(
      
      filename = function() {
        paste(input$treat_var, "csv", sep = ".")
      },
      content = function(file) {
        # Write to a file specified by the 'file' argument
        write.table(df_close(), file, sep = ",",row.names = FALSE)
      }
    )
    

# Total C -----------------------------------------------------------------

    # reading csv file
    # dat_totalc <- reactive({
    #   req(input$file_totalC)
    #   data <- read.csv(input$file_totalC$datapath)
    #   names <- unique(data$Treatment)
    #   # sizes <- sort(as.numeric(unique(data$Filter)))
    #   # num_parameters <- length(unique(data$Filter))
    #   
    #   list(data = data, names = names)
    # })
    
    
    # Choosing treatmentfg
    # observeEvent(input$source_totalC,{
    #   if(input$source_totalC==2){
    #     updateSelectInput(session, "treat_var_totalC", choices = input$treat1_totalC)
    #   }
    #   else{
    #     observeEvent(dat_totalc(), {
    #       updateSelectInput(session, "treat_var_totalC", choices = c(dat_totalc()$names,"ALL"))
    #       
    #     })
    #   }
    #   
    # })
    
    ## preparing data for regressions
    # raw.data.n_totalc<-reactive({
    #   req(input$file_totalC)
    #   raw.data = dat_totalc()$data
    #   raw.data = filter(raw.data, Distance>=0)
    #   # raw.data$Filter  = as.factor(raw.data$Filter)
    #   raw.data$Treatment  = as.factor(raw.data$Treatment)
    #   raw.data$Replicate  = as.factor(raw.data$Replicate)
    #   raw.data$Distance  = as.numeric(raw.data$Distance)
    #   
    #   ## normalizing mesocosm data with max value
    #   raw.data.n=raw.data %>% group_by(Treatment,Replicate) %>%
    #     mutate(max_c = max(Concentration,na.rm = TRUE))%>%
    #     mutate(norm_c = Concentration/max(Concentration,na.rm = TRUE)) %>%
    #     mutate(ln_c = ifelse((Concentration==0),NA,log(Concentration)))%>%
    #     mutate(ln_norm_c = ifelse((norm_c==0),NA,log(norm_c)))
    #   ## calculating the fractions based on normlaized values
    #   # raw.data.n = raw.data.n %>% group_by(Treatment,Replicate,Distance)%>%
    #   #   mutate(frac = norm_c/sum(norm_c)) %>%ungroup()
    #   
    #   raw.data.n
    # })
    
    # denorm_factor_totalc<- reactive({
    #   req(input$treat_var_totalC)
    #   raw.data.n_totalc()%>%ungroup()%>%
    #     subset(Treatment==input$treat_var_totalC)%>%
    #     summarise(average_maxC=mean(max_c,na.rm = TRUE))
    # })
    
    max_factor_totalc<- reactive({
      
      # if(input$data_type==1){
      req(input$file1)
      req(input$treat_var)
      raw.data.n()%>%ungroup()%>%
        subset(Treatment==input$treat_var)%>%
       summarise(Replicate=unique(Replicate),
                 maxC=unique(max_c,na.rm = TRUE))%>%
      mutate(lnmax=log(maxC))
      # }
      # else{
      #   req(input$treat_var)
      #   raw.data.n()%>%ungroup()%>%
      #     subset(Treatment==input$treat_var)%>%
      #     summarise(Replicate=unique(Replicate),
      #               maxC=unique(max_c,na.rm = TRUE))%>%
      #     mutate(lnmax=log(maxC))
      # }
        
    })
    
    ## Prediction forward 
    #single prediction table
    pred_totalC <-reactive({
      
      if(input$data_type==1){
      req(totalC())
      lambda =  table_var()$estimate
      lambdaup = c(table_var()$conf.low[1],table_var()$conf.high[2])
      lambdadown = c(table_var()$conf.high[1],table_var()$conf.low[2])
      # filt =  levels(table_pred()$Filter)
      # length(filt)
      tmax = 7*max(raw.data.n()$Distance) #input$time
      dt = 1
      time = seq(0,tmax,dt)
      cpred = (lambda[1] + time*lambda[2]) + mean(max_factor_totalc()$lnmax) #* as.numeric(denorm_factor_totalc())
      cpredup = (lambdaup[1] + time*lambdaup[2])+ max(max_factor_totalc()$lnmax) #* as.numeric(denorm_factor_totalc())
      cpreddown = (lambdadown[1] + time*lambdadown[2])+ min(max_factor_totalc()$lnmax) # * as.numeric(denorm_factor_totalc())
      # t = single_pred(tmax,dt,lambda,filters())
      t = tibble(time,cpred,cpredup,cpreddown)
      # t = cbind(time,cpred)
      }

      else{
        req(perSize())
        result_df <- data.frame()
        for (i in 1:length(filters())){
        # i=1
        lambda =  table_var()$estimate[(2*i-1):(2*i)]
        lambdaup = c(table_var()$conf.low[(2*i-1)],table_var()$conf.high[(2*i)])
        lambdadown = c(table_var()$conf.high[(2*i-1)],table_var()$conf.low[(2*i)])
        tmax = 7*max(raw.data.n()$Distance) #input$time
        dt = 1
        time = seq(0,tmax,dt)
        cpred = (lambda[1] + time*lambda[2]) + mean(max_factor_totalc()$lnmax) #* as.numeric(denorm_factor_totalc())
        cpredup = (lambdaup[1] + time*lambdaup[2]) + max(max_factor_totalc()$lnmax) #* as.numeric(denorm_factor_totalc())
        cpreddown = (lambdadown[1] + time*lambdadown[2]) + min(max_factor_totalc()$lnmax) # * as.numeric(denorm_factor_totalc())
        # t = single_pred(tmax,dt,lambda,filters())
        t = tibble(time,cpred,cpredup,cpreddown,Filter=filters()[i])
        result_df <- rbind(result_df,t)
        }
        result_df
      }
    })
    
    pred_distance_lod <-reactive({
      
      
      if(input$data_type==1){
      req(totalC())
      # browser()
      lambda =  table_var()$estimate
      lambdaup = c(table_var()$conf.low[1],table_var()$conf.high[2])
      lambdadown = c(table_var()$conf.high[1],table_var()$conf.low[2])
      
      # lod = log(input$LOD/as.numeric(denorm_factor_totalc())) # normalizing the LOD with an average factor across replicates
      lod = log(input$LOD) # normalizing the LOD with an average factor across replicates
      
      dpred = (lod - lambda[1] - mean(max_factor_totalc()$lnmax)) / lambda[2]
      dpredup = (lod - lambdaup[1] - max(max_factor_totalc()$lnmax)) / lambdaup[2]
      dpreddown = (lod - lambdadown[1] - min(max_factor_totalc()$lnmax)) / lambdadown[2]
      
      if(dpred<0){dpred=NaN}
      if(dpredup<0){dpredup=NaN}
      if(dpreddown<0){dpreddown=NaN}
      # t = single_pred(tmax,dt,lambda,filters())
      t = tibble(avg = round(dpred),upper = round(dpredup),lower = round(dpreddown))
      # t = cbind(time,cpred)
      }
      else{
        req(perSize())
        result_df <- data.frame()
        for (i in 1:length(filters())){
          # i=1
          lambda =  table_var()$estimate[(2*i-1):(2*i)]
          lambdaup = c(table_var()$conf.low[(2*i-1)],table_var()$conf.high[(2*i)])
          lambdadown = c(table_var()$conf.high[(2*i-1)],table_var()$conf.low[(2*i)])
          
          lod = log(input$LOD) # normalizing the LOD with an average factor across replicates
          
          dpred = (lod - lambda[1] - mean(max_factor_totalc()$lnmax)) / lambda[2]
          dpredup = (lod - lambdaup[1] - max(max_factor_totalc()$lnmax)) / lambdaup[2]
          dpreddown = (lod - lambdadown[1] - min(max_factor_totalc()$lnmax)) / lambdadown[2]
          
          if(dpred<0){dpred=NaN}
          if(dpredup<0){dpredup=NaN}
          if(dpreddown<0){dpreddown=NaN}
          # t = single_pred(tmax,dt,lambda,filters())
          t = tibble(Filter=filters()[i],avg = formatC(dpred,format="d"),
                     upper = formatC(dpredup,format="d"),lower = formatC(dpreddown,format="d"))
          result_df <- rbind(result_df,t)
        }
        result_df
      }
    })
    
    # plot of regressions for one selected treatment
    output$plot_pred_totalC <- renderPlot({
      # req(totalC())
      req(input$file1)
         
      
      aa = pred_totalC()%>%
        pivot_longer(cols = c(cpred,cpredup,cpreddown),values_to = "cpred", names_to = "cat")
        
        
        gg<- ggplot(aa)+
        geom_line(aes(x=time, y=cpred, color=cat),size=2.5,shape=21)+
        # geom_line(aes(x=time, y=cpredup),size=2.5,shape=21)+
        # geom_smooth(aes(x=Distance, y=ln_norm_c),method='lm')+
        # facet_grid(cols = vars(Filter))+
        theme_bw()+
        theme(text = element_text(size = 15))+
        ylab('ln(copies/ml)')+
        labs(x='distance (m)')+
        geom_hline(aes(yintercept = log(input$LOD), linetype = "LOD"),color = "red") +
        scale_linetype_manual(values = "dashed", name = " ")+
        scale_color_manual(values=c("black","blue","red"),name = " ",
                           breaks = c("cpred", "cpreddown", "cpredup"),
                           labels = c("avg", "lower", "upper"))+
        geom_point(data=data_regre(),aes(x=Distance, y=ln_c),size=2.5,shape=21)
      
      if (input$data_type==2) {
        gg <- gg + facet_grid(cols = vars(Filter))
      }
        
      if(input$data_in==2){
        gg<-gg + labs(x='Time (hr)')
      }
      
      gg
    }) 
    

    output$raw_totalc <- renderTable({
      
      # raw.data.n_totalc()%>%ungroup()%>%subset(Treatment==input$treat_var_totalC)%>%
      #   summarise(average_maxC=mean(max_c,na.rm = TRUE))
      # denorm_factor_totalc()
      
      pred_distance_lod()
      
      # as.numeric(denorm_factor_totalc())
      # raw.data.n()%>%subset(Treatment==input$treat_var)%>%
      #   mutate(average_maxC=mean(max_c,na.rm = TRUE))##%>%
      # subset(Filter==10)
      # pred_totalC()
      
    },digits=1)
    
    output$pred_totalc <- renderTable({
      # pred_totalC()
      max_factor_totalc()
    })
    
# SMIM --------------------------------------------------------------------
    
    observeEvent(input$file1,{
      if(input$data_type==2){
        updateSelectInput(session, "size2smim", choices = filters())
      }
    })
    
    
    
    logT1 <- 0
    logT2 <- 5
    
    C0 <- reactive({
      if(input$data_type==2){
        data = data_regre()%>%subset(Filter==input$size2smim)%>%
          subset(Distance==0)
        C0 = mean(data$Concentration)
      }
      else{
       C0 = mean(max_factor_totalc()$maxC)
      }
    })
    
    output$SMIM_Plot <- renderPlot({
      C0 = C0() #mean(max_factor_totalc()$maxC) #input$C0
      rs = input$rss #optim_param()$par[1] #input$rs
      rh = optim_param()$par[1] #input$rh
      U = input$Us; D = input$Ds; lambda=input$lambdas
      beta=input$betas
      # U = optim_param()$par[3] #input$U; 
      # D = optim_param()$par[4] #input$D; 
      # lambda=optim_param()$par[5] #input$lambda
      # beta=optim_param()$par[6] #input$beta;
      # logT1=optim_param()$par[7] #input$logT1;
      # logT2=optim_param()$par[8] #input$logT2
      xlim = input$xlim
      
      C0s = C0() # mean(max_factor_totalc()$maxC) #input$C0
      rss = input$rss
      rhs = input$rhs
      Us= input$Us; Ds = input$Ds; lambdas=input$lambdas
      betas=input$betas;logT1s=logT1;logT2s=logT2
      xlims = input$xlim
      
      
      if(input$data_type==2){
        
      data= data_regre()%>%subset(Filter==input$size2smim) 
      }
      else{
        data = data_regre()
      }
      
      cols <- c('Optim' = "black", 'Man' = "red")
      dfbase = C_func(xlim,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
      dfsens = C_func(xlims,C0s,rss,rhs,Us,Ds,lambdas, betas,logT1s,logT2s)
      ggplot(dfbase)+geom_line(aes(x,C,colour='Optim'),size=2)+
        geom_line(data=dfsens,aes(x,C,color='Man'),linetype = "dashed",size=2)+
        theme_bw()+theme(text = element_text(size = 15))+
        scale_y_continuous(name='Concentration(copies/ml)')+
        scale_color_manual(name = "Model",values=cols)+
        geom_point(data=data,aes(x=Distance, y=Concentration),size=2.5,shape=21)+
        geom_vline(aes(xintercept =input$xloc), linetype = "dashed",color = "blue")+
        annotate("text", x = input$xloc, y = C0, label = textplot(), size = 4, color = "black")+
        annotate("text", x = input$xloc, y = .95*C0, label = textplot2(), size = 4, color = "red")
      #norm_c*input$C0
    })
    
    
    textplot <- reactive({
      C0 = mean(max_factor_totalc()$maxC) #input$C0
      rs = input$rss # optim_param()$par[1] #input$rs
      rh = optim_param()$par[1] #input$rh
      U = input$Us; D = input$Ds; lambda=input$lambdas
      beta=input$betas
      # U = optim_param()$par[3] #input$U; 
      # D = optim_param()$par[4] #input$D; 
      # lambda=optim_param()$par[5] #input$lambda
      # beta=optim_param()$par[6] #input$beta;
      # logT1=optim_param()$par[7] #input$logT1;
      # logT2=optim_param()$par[8] #input$logT2
      xloc = input$xloc
      Cloc = C_loc(xloc,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
      
      paste("Optim case: ",
            "C =",round(Cloc,2),"copies/ml")
    })
    
    textplot2 <- reactive({
      C0 = mean(max_factor_totalc()$maxC) #input$C0
      rs = input$rss
      rh = input$rhs
      U = input$Us; D = input$Ds; lambda=input$lambdas
      beta=input$betas #;logT1=input$logT1s;logT2=input$logT2s
      xloc = input$xloc
      Cloc = C_loc(xloc,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
      paste("Man case: ",
            "C =",round(Cloc,2),"copies/ml")
    })
    
    # output$text_smim_2 <- renderText({
    #   C0 = mean(max_factor_totalc()$maxC) #input$C0
    #   rs = input$rss
    #   rh = input$rhs
    #   U = input$Us; D = input$Ds; lambda=input$lambdas
    #   beta=input$betas;logT1=input$logT1s;logT2=input$logT2s
    #   xloc = input$xloc
    #   Cloc = C_loc(xloc,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
    #   paste("Man case:\n",
    #         "At x=",input$xloc,"m :",
    #         "C =",round(Cloc,2),"copies/ml")
    # })
    
    output$koptim <- renderText({
      # C0 = input$C0
      rs = input$rss # optim_param()$par[1] #input$rs
      rh = optim_param()$par[1] #input$rh
      U = input$Us; D = input$Ds; lambda=input$lambdas
      beta=input$betas
      # U = optim_param()$par[3] #input$U; 
      # D = optim_param()$par[4] #input$D; 
      # lambda=optim_param()$par[5] #input$lambda
      # beta=optim_param()$par[6] #input$beta;
      # logT1=optim_param()$par[7] #input$logT1;
      # logT2=optim_param()$par[8] #input$logT2
      
      k = k_func(rs,rh,U,D,lambda, beta,logT1,logT2)
      paste('<b>',round(k,4),'</b>')
    })
    
    output$ks <- renderText({
      # C0 = input$C0
      rs = input$rss
      rh = input$rhs
      U = input$Us; D = input$Ds; lambda=input$lambdas
      beta=input$betas#;logT1=input$logT1s;logT2=input$logT2s
      k = k_func(rs,rh,U,D,lambda, beta,logT1,logT2)
      paste('<b>',round(k,4),'</b>')
    })
    
    output$img <- renderImage({
      
      list(src = "SMIM_schematic.png",
           width = "75%",
           height = 300)
      
    }, deleteFile = F)
#finding parameters(optimization)
    
    optim_param <-reactive({
      req(input$file1)
      #rs,rh,U,D,lambda, beta,logT1,logT2
      par = c(1e-3,1e-2)
      
      fr<- function(x) {
        
        if(input$data_type==2){
          data= table_var()%>%subset(Filter==input$size2smim)
          kdata = -data$estimate[2]
        }
        else{
          kdata = -table_var()$estimate[2]
        }
        
        
        
        U = input$Us; D = input$Ds; lambda=input$lambdas
        beta=input$betas
        rs = input$rss
        abs(k_func(rs,x[1],U,D,lambda, beta,logT1,logT2)-kdata)
        # abs(k_func(x,1e-3,1,3e-2,0.1,0.9,0,5)-kdata)
        
      }
      # kdata = table_var_totalc()$estimate[2]
      xoptim = optim(par,fr,lower = rep(0, 1), method = "SANN")
      
    })
    
    output$whatever_optim<- renderTable({
      kdata = -table_var()$estimate[2]
      c(kdata,optim_param()$value)
    },digits = 3)
    
    # output$rs <- renderText({
    #   paste('<b>',round(optim_param()$par[1],3),'</b>')
    # })
   
    output$rh <- renderText({
      paste('<b>',formatC(optim_param()$par[1],format="e",digits=1),'</b>')
    })
    # 
    # output$U <- renderText({
    #   paste('<b>',round(optim_param()$par[3],3),'</b>')
    # })
    # 
    # output$D <- renderText({
    #   paste('<b>',round(optim_param()$par[4],3),'</b>')
    # })
    # 
    # output$lambda <- renderText({
    #   paste('<b>',round(optim_param()$par[5],3),'</b>')
    # })
    # 
    # output$beta <- renderText({
    #   paste('<b>',round(optim_param()$par[6],3),'</b>')
    # })
    
    # output$logT1 <- renderText({
    #   paste('<b>',round(optim_param()$par[7],3),'</b>')
    # })
    # 
    # output$logT2 <- renderText({
    #   paste('<b>',round(optim_param()$par[8],3),'</b>')
    # })
    
    


           
} 

# end ---------------------------------------------------------------------


# Run the app
shinyApp(ui, server)
