library(shiny)
library(ggplot2)
library(bslib)
library(tidyverse)
library(glue)
library(colourpicker)
library(patchwork)
library(gplots)
library(pheatmap)
library(ggbeeswarm)
library(DT)
library('fgsea')
library('GSEABase')


draw_table_Summary<-function(metadata){
  
  col1<-c('Age', 'Dianosis', 'Postmortem interval(hours)', 'onset-age', 'Disease duration(years)')
  col2<-c('double', 'factor', 'double', 'double', 'double')
  row_age<-paste(round(mean(pull(meta_data,'age_of_death')),2),"(+/-",round(sd(pull(meta_data,'age_of_death')),2),")")
  row_Dig<-"Normal, Huntington'Disease"
  row_PMI<-paste(round(mean(pull(meta_data,'pmi'),na.rm=TRUE),2),"(+/-",round(sd(pull(meta_data,'pmi'),na.rm=TRUE),2),")")
  row_onset<-paste(round(mean(pull(meta_data,'age_of_onset'),na.rm=TRUE),2),"(+/-",round(sd(pull(meta_data,'age_of_onset'),na.rm=TRUE),2),")")
  row_dd<-paste(round(mean(pull(meta_data,'Duration'),na.rm=TRUE),2),"(+/-",round(sd(pull(meta_data,'Duration'),na.rm=TRUE),2),")")
  col3<-c(row_age, row_Dig, row_PMI, row_onset, row_dd)
  #generate the tibble for UI
  tbl <- tibble(
    'Column Name' = col1,
    'Type'=col2,
    'Mean(sd) or distinct values'=col3
  )
  #print(tbl)
  return(tbl)
}


#this funtion is used to sort the metadata table
draw_table_metadata<-function(metadata){
  #print(metadata)
  data<-metadata %>% dplyr::select(
    Run,
    age_of_death,
    Diagnosis,
    mrna.seq_reads,
    pmi,
    rin
  ) 
  table <- DT::datatable(data, extensions = 'Buttons', class = "display",
                         options = list(paging = TRUE, searching = TRUE,
                                        fixedColumns = TRUE, autoWidth = TRUE,
                                        ordering = TRUE, dom = 'Bfrtip',
                                        buttons = c('copy', 'csv')))
  return(table)
}

#This function is used to draw the violin plot
draw_violin_metadata<-function(metadata, y_axis){
  plot<- ggplot(metadata) +
          geom_violin(aes(x=Diagnosis,y=get(y_axis),fill=Diagnosis),na.rm=TRUE) +
          ylab(y_axis)
  return(plot)
}


#function to judge the point is passed?
ifpassed <- function(p1,p2){
  return(p1 && p2)
}


#this function is used to draw the summary table of counts matrix
draw_summary_matrix <- function(matrix){
  gene_num=nrow(matrix)
  sample_num=length(matrix)-7
  pass_num=sum(matrix['passed']==TRUE)
  pass_per=round(pass_num/gene_num*100,2)
  nonpass_num=gene_num-pass_num
  nonpass_per=round(nonpass_num/gene_num*100,2)
  col1=c("The number of Genes", "The number of Samples", "The number of gene passing filter", "The number of gene not passing filter")
  col2=c(gene_num, sample_num, paste(pass_num,"(", pass_per, "%)"), paste(nonpass_num,"(", nonpass_per, "%)"))
  tbl<-tibble(
    Category=col1,
    Number_percentage=col2
  )
  return(tbl)
}


#this function is used to draw the scatter plot of counts
draw_scatter_matrix <- function(matrix){
  plot1<-ggplot(matrix,aes(x=mid, y=var, col=passed))+
          geom_point()+
          ggplot2::scale_y_log10()+
          ggplot2::scale_x_log10(oob = scales::squish_infinite)
  plot2<-ggplot(matrix,aes(x=mid, y=non_zero, col=passed))+
          geom_point()+
          ggplot2::scale_x_log10(oob = scales::squish_infinite)
  plot<-plot1+plot2
  return(plot)
}


#this function is used to draw the heatmap for counts matrix
draw_heatmap_matrix <- function(matrix){
  data<- filter(matrix, 
            passed==TRUE
         )
  n_col=ncol(data)
  data <- data[,2:(n_col-6)] %>% as.matrix()
  #print(head(data))
  #print(meta_data['Diagnosis'])
  color<-if_else(meta_data['Diagnosis']=="Huntington's Disease", 'red', 'blue')
  #print(color)
  plot<-heatmap.2(
    data,
    ColSideColors=color,
    scale='row',
    symkey=FALSE, density.info="none", trace="none",
    ylab="Genes", 
    xlab="Samples",
    main="The heatmap for count matrix",
    labRow=FALSE
  )
  return(plot)
}

#this function is used to generate the PCA projection
draw_PCA_matrix <- function(matrix,meta,PC1, PC2, topPC){
  #print(PC1)
  #print(PC2)
  data<- filter(matrix, 
                passed==TRUE
          )
  n_col=ncol(data)
  data_PCA <- data[,2:(n_col-6)] %>% as.matrix()
  rownames(data_PCA) <- data$X
  pca <- prcomp(
    t(data_PCA),
    center=FALSE, # prcomp centers features to have mean zero by default, but we already did it
    scale=TRUE # prcomp scales each feature to have unit variance
  )
  
  #print(pca$x[,1:topPC])
  #plot the first beeswarm 
  plot1<- as_tibble(pca$x[,1:topPC]) %>%
    pivot_longer(everything(),names_to="PC",values_to="projection") %>%
    mutate(PC=fct_relevel(PC,str_c("PC",1:20))) %>%
    ggplot(aes(x=PC,y=projection)) +
    geom_beeswarm() + labs(title="PCA Projection Plot") +
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
  
  
  
  #plot the second plot
  plot2_data<- tibble(
    type=meta$Diagnosis,
    PC1=pca$x[,PC1],
    PC2=pca$x[,PC2]
  )
  #print(plot2_data)
  title="Counts PCA"
  var_p <- pca$sdev^2 / sum( pca$sdev^2 )
  x_title <- paste0("PC",PC1," :", round(var_p[PC1] * 100),"% variance")
  y_title <- paste0("PC",PC2," : ", round(var_p[PC2] * 100), "% variance")
  plot2<- as_tibble(plot2_data) %>%
          ggplot(aes(x=PC1,y=PC2, col=type)) +
          geom_point()+
          ggtitle(title)+
          xlab(x_title)+
          ylab(y_title)
  
  return(plot1/plot2)
}


#this function is used to return the DE table
draw_table_DE<- function(dataf, slider) {
  data<-dataf %>%
    filter(
      log10(padj)<slider
    ) %>%
    mutate(
      pvalue=formatC(pvalue,format = 'e'),
      padj=formatC(padj,format = 'e')
    )
  
  table <- DT::datatable(data, extensions = 'Buttons', class = "display",
                         options = list(paging = TRUE, searching = TRUE,
                                        fixedColumns = TRUE, autoWidth = TRUE,
                                        ordering = TRUE, dom = 'Bfrtip',
                                        buttons = c('copy', 'csv')))
  return(table)
}

#this function is for drwaing the volcano plot
volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    data<- dataf %>%
      filter(
        !is.na(padj)
      )%>%
      mutate(
        'colorf'=ifelse(log10(padj)<slider,TRUE, FALSE)
        
      ) 
    plot <- ggplot(data,aes( !!sym(x_name), -log10(!!sym(y_name)) , color=colorf) )+
      geom_point() +
      scale_color_manual(values=c(color1,color2))+
      labs(color=glue("pvalue<10^{slider}")) 
    
    return(plot)
  }


#this function is used to draw the topest gsea result
draw_top_gsea <- function(fgsea_result, num){
  fgsea_result %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES))
    
  plot <- ggplot(fgsea_results[1:num,]) +
    geom_bar(aes(x=pathway, y=NES, fill = NES > 0), stat='identity') +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
    theme_minimal() +
    ggtitle('fgsea results for Hallmark MSigDB gene sets') +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') +
    coord_flip()
  return(plot)
}

#this function is used to return the table of gsea 
draw_table_gsea <-function(fgsea_result, pval_gsea,thre_gsea){
  #("Positive", "Negative", "All")
  if (thre_gsea=="Positive"){
    up<-10000
    down<-0
  }
  if (thre_gsea=="Negative"){
    up<-0
    down<-(-10000)
  }
  if (thre_gsea=="All"){
    up<-10000
    down<-(-10000)
  }
  
  data<-fgsea_result %>% 
        filter(
          log10(padj)<pval_gsea,
          NES>down,
          NES<up
        )
  
  table <- DT::datatable(data, extensions = 'Buttons', class = "display",
                         options = list(paging = TRUE, searching = TRUE,
                                        fixedColumns = TRUE, autoWidth = TRUE,
                                        ordering = TRUE, dom = 'Bfrtip',
                                        buttons = c('copy', 'csv')))
  return(table)
}

#this function is used to draw the second plot for gsea
draw_plot2_gsea <- function(fgsea_result, pval_plot_gsea){
  data<-fgsea_result %>% 
    mutate(
      passed=log10(padj)<pval_plot_gsea,
    )
  
  #draw the plot
  plot<-ggplot(data,aes(x=NES , y=padj, col=passed))+
    geom_point()+
    ggplot2::scale_y_log10()
  return(plot)
}

ui <- fluidPage(
  #head title
  titlePanel("BF591 Final Project"),
  p("Data set: GSE64810"),
  p("publication: The caudate nucleus undergoes dramatic and unique transcriptional changes in human prodromal Huntington's disease brain. BMC Med Genomics 2019 Oct 16;12(1):137. PMID: 31619230"),
  #the main four tabs
  tabsetPanel(
    #tab1 samples
    tabPanel(
      "Samples",
      sidebarLayout(
        sidebarPanel(
          style="height:600px;",
          fileInput("Sample_file", paste0("Sample file"), accept='.csv'),
          actionButton("Submit_1",label='Submit', style="background-color: #9df2fa;",width=180)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Summary",
              tableOutput("Summary")
            ),
            tabPanel(
              "Table",
              dataTableOutput("metadata")

            ),
            tabPanel(
              "Plots",
              sidebarLayout(
                sidebarPanel(
                  radioButtons('meta_y_axis',"Y-axis is",
                               choices = c("age_of_death", "mrna.seq_reads","pmi","rin"), selected="age_of_death"),
                  p("pmi=Postmortem interval"),
                  width=3
                ),
                mainPanel(
                  plotOutput("Summary_plot"),
                  width=8
                )
              )
            ),
          )
        )
      )
    ),
    #tab2 counts
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(
          style="height:600px;",
          fileInput("matrix", paste0("Normalized counts matrix"), accept='.csv'),
          sliderInput(inputId = "var_per", min = 0, max = 1, 
                      label = "The percentile of variance for genes included", value = 0.95, step = 0.05),
          sliderInput(inputId = "non_zero", min = 0, max = 69, 
                      label = "The least number of non-zero samples for genes included", value = 69, step = 1),
          actionButton("Submit_2",label='Submit', style="background-color: #9df2fa;",width=180)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Summary",
              tableOutput("Summary_Counts")
            ),
            tabPanel(
              "scatter_plots",
              plotOutput("scatter_Counts")
            ),
            tabPanel(
              "heatmap",
              plotOutput("heatmap_Counts")
            ),
            tabPanel(
              "PCA_plot",
              sidebarLayout(
                sidebarPanel(
                  sliderInput(inputId = "TopN", min = 1, max = 20, 
                              label = "the top N principal components for beeswarm plot", value = 5, step = 1),
                  sliderInput(inputId = "PC1", min = 1, max = 20, 
                              label = "The first principal components you want to plot", value = 1, step = 1),
                  sliderInput(inputId = "PC2", min = 1, max = 20, 
                              label = "The second principal components you want to plot", value = 2, step = 1),
                  width=4
                ),
                mainPanel(
                  style="height:900px;",
                  plotOutput("PCA_Counts"),
                  width=8
                )
              )
            )
          )
        )
      )
    ),
    #tab3 differential expr
    tabPanel(
      "DE",
      sidebarLayout(
        sidebarPanel(
          fileInput("deseq", paste0("Load differential expression results"), accept='.csv'),
          radioButtons('xaxis',"Choose the column for the x-axis",
                       choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"), selected="log2FoldChange"),
          radioButtons('yaxis',"Choose the column for the y-axis",
                       choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"), selected="padj"),
          colourInput('base', "Base point color", showColour="both", value="#22577A"),
          colourInput('highlight', "Highlight point color", showColour="both", value="#FFCF56"),
          sliderInput(inputId = "slider", min = -30, max = 0, 
                      label = "Select the magnitude of the p adjusted coloring:", value = -10, step = 1),
          actionButton("plot_deseq",label='Refresh', icon = icon("refresh","fa"), style="background-color: #9df2fa;",width=180)
          
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Table",
                     dataTableOutput("DE_table")  #need a search tool and the table
            ),
            tabPanel("Plot",
                     plotOutput("volcano")
            )
          )
        )
      )
    ),
    #tab4 GSEA
    tabPanel(
      "GSEA",
      sidebarLayout(
        sidebarPanel(
          p("Gene Set Enrichment Analysis"),
          fileInput("deseq_GSEA", paste0("Load differential expression results"), accept='.csv'),
          radioButtons('metric_GSEA',"Choose the ranking metric for the GSEA",
                       choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"), selected="log2FoldChange"),
          #the parameter for fgsea, minSize and maxSize
          sliderInput(inputId = "minsize", min = 0, max = 50, 
                      label = "Select the parameter(minSize) for fgsea", value = 15, step = 1),
          sliderInput(inputId = "maxsize", min = 400, max = 600, 
                      label = "Select the magnitude of the p adjusted coloring:", value = 500, step = 5),
          actionButton("submit_GSEA",label='Upload/update', icon = icon("upload", "fa"),width=180),
          actionButton("download_GSEA",label='Save the result', icon = icon("download", "fa"),width=180)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Top Results",
                     sidebarLayout(
                       sidebarPanel(
                         style="height:600px;",
                         sliderInput(inputId = "Top_num_gsea", min = 1, max = 50, 
                                     label = "Select the number of top pathways for ploting", value = 5, step = 1),
                         actionButton("Top_gsea_update",label='Refresh', icon = icon("refresh","fa"), style="background-color: #9df2fa;",width=180)
                       ),
                       mainPanel(
                         plotOutput("top_plot")
                       )
                     )
            ),
            tabPanel("Table",
                     sidebarLayout(
                       sidebarPanel(
                         style="height:600px;",
                         sliderInput(inputId = "pval_gsea", min = -30, max = 0, 
                                     label = "Select the magnitude of the p adjusted coloring:", value = -10, step = 1),
                         radioButtons(inputId='thre_gsea',"Choose the NES score(+/-): ",
                                      choices = c("Positive", "Negative", "All"), selected="All"),
                         actionButton("table_gsea_update",label='Refresh', icon = icon("refresh","fa"), style="background-color: #9df2fa;",width=180)
                       ),
                       mainPanel(
                         dataTableOutput("GSEA_table")  
                       )
                     )
            ),
            tabPanel("Plot",
                     sidebarLayout(
                       sidebarPanel(
                         style="height:600px;",
                         sliderInput(inputId = "pval_plot_gsea", min = -30, max = 0, 
                                     label = "Select the magnitude of the p adjusted coloring:", value = -10, step = 1),
                         actionButton("plot_gsea_update",label='Refresh', icon = icon("refresh","fa"), style="background-color: #9df2fa;",width=180)
                       ),
                       mainPanel(
                         plotOutput("plot_gsea")
                       )
                     )
            )
          )
        )
      )
    ),
  )
)
#input file size setting
options(shiny.maxRequestSize=30*1024^2)

server <- function(input, output, session) {
  #load data 1, the meta data
  observeEvent(input$Submit_1,{
    if (length(input$Sample_file$datapath!=0)){
      meta_data<<-read.csv(input$Sample_file$datapath, header = TRUE) %>%
        as.tibble()
      #print(summary(meta_data))
      return(meta_data)
    }
  })
  #listen to the second button 
  observeEvent(input$Submit_2,{
    if (length(input$matrix$datapath!=0)){
      matrix<-read.csv(input$matrix$datapath, header = TRUE) %>%
              as_tibble()
      #print(nrow(matrix))
      len_matrix=length(matrix)
      matrix_filtered<<- matrix%>%
                mutate(
                  var=apply(matrix[,2:len_matrix], 1, var),
                  mid=apply(matrix[,2:len_matrix], 1, median),
                  non_zero=rowSums(matrix != 0)-1
                )%>%
                mutate(
                  passed_1=(var>=quantile(var, input$var_per)),
                  passed_2=(non_zero>=input$non_zero),
                  passed=mapply(ifpassed,passed_1,passed_2)
                )
      #print(head(matrix_filtered))
      return(matrix_filtered)
    }
  })
  #listen to the third button
  observeEvent(input$plot_deseq,{
    if (length(input$deseq$datapath!=0)){
      DE_data<<-read.csv(input$deseq$datapath, header = TRUE) %>%
        as.tibble()
      return(DE_data)
    }
  })
  #listen to the GSEA submit buttion
  observeEvent(input$submit_GSEA,{
    if (length(input$deseq_GSEA$datapath!=0)){
      GSEA_data<<-read.csv(input$deseq_GSEA$datapath, header = TRUE) 
    }
    GSEA_data<<-GSEA_data%>% as.tibble() %>%
      arrange(input$metric_GSEA)
    #generate the ranking list
    rnks<-pull(GSEA_data,input$metric_GSEA)
    rnk_list <- setNames(rnks, pull(GSEA_data,'symbol'))
    #print("update data")
    #do GSEAnalysis
    hallmarks_gmt <- getGmt(con='h.all.v7.0.symbols.gmt')
    hallmark_pathways_GSEABase <- geneIds(hallmarks_gmt)
    fgsea_results <<- fgsea(hallmark_pathways_GSEABase, rnk_list, minSize = input$minsize, maxSize=input$maxsize) %>% 
                      as_tibble()
    print(fgsea_results)
    return(GSEA_data)
  })
  
  output$Summary <- renderTable({
    input$Submit_1
    return(draw_table_Summary(meta_data))
  })
  output$metadata <- DT::renderDataTable({
    input$refresh_meta
    return(draw_table_metadata(meta_data))
  })
  output$Summary_plot <- renderPlot({
    input$Submit_1
    return(draw_violin_metadata(meta_data, input$meta_y_axis))
  })
  output$Summary_Counts <-renderTable({
    input$Submit_2
    return(draw_summary_matrix(matrix_filtered))
  })
  output$scatter_Counts <- renderPlot({
    input$Submit_2
    return(draw_scatter_matrix(matrix_filtered))
  })
  output$heatmap_Counts <- renderPlot({
    input$Submit_2
    plot<-draw_heatmap_matrix(matrix_filtered)
    return(plot)
  })
  output$PCA_Counts <- renderPlot({
    input$Submit_2
    plot<-draw_PCA_matrix(matrix_filtered, meta_data, input$PC1, input$PC2, input$TopN)
    return(plot)
  },
  height=900
  )
  #the Differetial expression table
  output$DE_table<-DT::renderDataTable({
    input$plot_deseq
    return(draw_table_DE(DE_data,input$slider))
  })
  output$volcano <- renderPlot({
    volcano_plot(DE_data, input$xaxis, input$yaxis, input$slider, input$base, input$highlight)
  })
  #gsea top plot
  output$top_plot <- renderPlot({
    input$Top_gsea_update
    return(draw_top_gsea(fgsea_results,input$Top_num_gsea))
  })
  
  #gsea table
  output$GSEA_table<-DT::renderDataTable({
    input$table_gsea_update
    return(draw_table_gsea(fgsea_results, input$pval_gsea,input$thre_gsea))
  })
  
  #gsea plot 2
  output$plot_gsea <- renderPlot({
    input$plot_gsea_update
    return(draw_plot2_gsea(fgsea_results, input$pval_plot_gsea))
  })
  
}

shinyApp(ui = ui, server = server)